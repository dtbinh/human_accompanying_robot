# utility files
# originally from MPC
# modified for single target search in early 2016
# later modified for human companion robot 
# Chang Liu 6/7/16
# restart working on 6/28/17

module CarMpc

using JuMP
using Ipopt
using CarMpcUtils
using Polynomials, Interpolations

import JuMP.JuMPArray, JuMP.Variable, JuMP.NonlinearParameter

export MpcProblem, initializeMPCProblem, updateSolveMpcProblem
# export localizeVehicle, generateReference

type MpcProblem
  # Vehicle
  robot::Robot
  # obstacle
  human::Array{Obstacle} 
  obs::Array{Obstacle}
  # field
  field::Field
  # Constraints
  uLB::Array{Float64,1}
  uUB::Array{Float64,1}
  zLB::Array{Float64,1}
  zUB::Array{Float64,1}
  # Parameters
  nz::Int64
  nu::Int64
  tuning::Tuning
  # JuMP model variables
  model::Model
  z::Array{Variable,2}
  u::Array{Variable,2}
  # JuMP model parameters
  z0::Array{NonlinearParameter,1}
  # u0::Array{Float64,1}
  z_h::Array{NonlinearParameter,2}
end

### Set up MPC optimization problem
function initializeMPCProblem(robot::Robot, human::Array{Obstacle}, obs::Array{Obstacle}, field::Field, tuning::Tuning)
  # human states are from ROS message. obs are obstacle in the environment based on prior knowledge

  ### Parameters

  N = tuning.N
  dt = tuning.dt
  safe_dis = tuning.safe_dis
  safe_margin = tuning.safe_margin
  cmft_dis = tuning.cmft_dis

  uLB = [robot.steerMin, robot.accMin]
  uUB = [robot.steerMax, robot.accMax]
  zLB = [0.0;0.0]
  zUB = [10.0;10.0]
  maxV = robot.maxV
  init_pose = robot.z

  # determine the number of humans and size of z_h
  numH = size(human,1); 
  if (numH == 0)
    zhSize = 0;
  else
    zhSize = numH*length(human[1].z);
  end

  # determine the number of static obstacles
  numObs = size(obs,1);

  ### Optimization problem
#   model = Model(solver=IpoptSolver(linear_solver="ma27",
#                                      max_iter=100,
#                                      max_cpu_time=dtMPC))
  model = Model(solver=IpoptSolver(max_iter=100, max_cpu_time=100.0, print_timing_statistics="no")) #max_cpu_time=dtMPC, hessian_approximation ="limited-memory"

  ### Dimensions
  nz = 4
  nu = 2

  ### Dummies
  # Parameters for optimization problem that can be modified online
  @NLparameter(model, z0[i=1:nz] == init_pose[i]) #z0 = init_pose
  # @NLparameter(model, u0 == zeros(nu)) #u0 = zeros(nu)

  @NLparameter(model, z_h[i=1:zhSize,j=1:N+1] == 0) #z_h = zeros(zhSize,N+1)

  ### Variables
  @variable(model, z[1:nz,1:N+1])
  @variable(model, u[1:nu,1:N])
  @variable(model, dummy_t[1:N+1]>=0)

  # ### initial solution
  # ux_init = (ter_pos[1]-z0[1])/N;
  # uy_init = (ter_pos[2]-z0[2])/N;

  # if (ux_init^2+uy_init^2) <= maxV^2
  #   U_init = [ux_init*ones(1,N);uy_init*ones(1,N)]
  # else
  #   angle = atan2(ter_pos[2]-z0[2],ter_pos[1]-z0[1])
  #   U_init = [maxV*cos(angle)*ones(1,N);maxV*sin(angle)*ones(1,N)]
  # end

  # Z_init = zeros(2,N)
  # for k=1:N
  #   if k == 1
  #     Z_init[1,k] = z0[1] + dt*U_init[1,k]
  #     Z_init[2,k] = z0[2] + dt*U_init[2,k]
  #     # println(Z_init)
  #   else
  #     Z_init[1,k] = Z_init[1,k-1] + dt*U_init[1,k]
  #     Z_init[2,k] = Z_init[2,k-1] + dt*U_init[2,k]
  #     # println(Z_init)
  #   end
  # end

  # for k = 1:N
  #   setValue(u[k],U_init[k])
  #   setValue(z[k],Z_init[k])
  # end

  ### Objective function

  # @defNLExpr(obj[j=1:N], abs((z[1,j]-ter_pos[1])^2+(z[2,j]-ter_pos[2])^2-ter_r^2))
  # @setNLObjective(model, Min, sum{obj,i=1:N})

  # @setNLObjective(model, Min, sum{abs((z[1,j]-ter_pos[1])^2+(z[2,j]-ter_pos[2])^2-ter_r^2),j=1:N})
  @NLobjective(model, Min, sum(dummy_t[i] for i in 1:N)) # abs((z[1,k]-z_h[1])^2+(z[2,k]-z_h[2])^2-cmft_dis^2) + 2*(z[4,k]-v_h)^2)

  ### Constraints
  # initial condition
  @NLconstraint(model, initconstr[i=1:nz], z[i,1] == z0[i]) # initconstr is the name of constraints z[i.1]==z0[i]

  for k=1:N+1
    ### Dynamics
    if k <= N
    @NLconstraints(model, begin
      z[1,k+1] == z[1,k] + dt*z[4,k]*cos(z[3,k])
      z[2,k+1] == z[2,k] + dt*z[4,k]*sin(z[3,k])
      z[3,k+1] == z[3,k] + dt*u[1,k]
      z[4,k+1] == z[4,k] + dt*u[2,k]
      end)
    end

    ### state and input constraints
    @constraint(model, 0 <= z[4,k] <= maxV)
    # stay within the field boundary
    # @constraint(model, zLB[1] <=z[1,k] <= zUB[1])
    # @constraint(model, zLB[2] <=z[2,k] <= zUB[2])  

    if k <= N
      @constraint(model, uLB[1] <=u[1,k] <= uUB[1])
      @constraint(model, uLB[2] <=u[2,k] <= uUB[2])
    end

    ### epigraph variable for obj
    @NLconstraint(model,dummy_t[k]>=(z[1,k]-z_h[1,k])^2 + (z[2,k]-z_h[2,k])^2 - cmft_dis^2 + 2*(z[4,k]-z_h[3,k])^2)
    @NLconstraint(model,dummy_t[k]>=-((z[1,k]-z_h[1,k])^2 + (z[2,k]-z_h[2,k])^2 - cmft_dis^2) + 2*(z[4,k]-z_h[3,k])^2)

    ### collision avoidance    
    # avoid static obstacles
    for j = 1:numObs
      if obs[j].shape == "circle" # circle
        @NLconstraint(model,(obs[j].z[1]-z[1,k])^2+(obs[j].z[2]-z[2,k])^2 >= obs[j].size^2);
      end
      if obs[j].shape == "ellips" # ellipsoid
        # to fill
      end
    end

    # avoid human
    @NLconstraint(model,(z[1,k] - z_h[1,k])^2 + (z[2,k] - z_h[2,k])^2 >= safe_dis^2)
  end

  ### Solve dummy problem
  solve(model)

  return MpcProblem(robot, human, obs, field, uLB, uUB, zLB, zUB, nz, nu, tuning, model, z, u, z0, z_h)

end

### Update and solve MPC optimization problem
# function updateSolveMpcProblem(mpc::MpcProblem, ZRef::Array{Float64,2},
#                                URef::Array{Float64,2}, z0::Array{Float64,1},
#                                u0::Array{Float64,1})
function updateSolveMpcProblem(mpc::MpcProblem, z0::Array{Float64,1}, z_h::Array{Float64,2})
  # Parameters
  nz, nu, dt, N = mpc.nz, mpc.nu, mpc.tuning.dt, mpc.tuning.N

  # Warm start
  UPred = hcat(getvalue(mpc.u[:,2:N]), getvalue(mpc.u[:,N]))
  ZPred = simRobotModel(mpc.robot, UPred, dt)
  # map(setvalue, collect(mpc.z[:,:]), collect(ZPred[:,:]))
  # map(setvalue, collect(mpc.u[:,:]), collect(UPred[:,:]))
  setvalue(mpc.z, ZPred)
  setvalue(mpc.u, UPred)

  # Update model parameters
  # mpc.z0[:] = z0
  # mpc.z_h[:] = z_h
  # map(setvalue, collect(mpc.z0[:]), collect(z0[:]))
  # map(setvalue, collect(mpc.z_h[:,:]), collect(z_h[:,:]))
  setvalue(mpc.z0, z0)
  setvalue(mpc.z_h, z_h)

  # Solve updated problem
  tStart = time()
  solve(mpc.model)
  solveTime = time() - tStart

  return getvalue(mpc.z), getvalue(mpc.u), solveTime
end

### Localize vehicle in current lane on field and compute reference speed
# function localizeVehicle(zCurr::Array{Float64,1},field::TrackMap)
#   # Parameters
#   nodes, wayPointers, nLanes = field.nodes, field.wayPointers, field.nWays
#   speed = field.speed
#   nNodesPolyFront, nNodesPolyBack, nNodesThreshold = field.nNodesPolyFront,
#     field.nNodesPolyBack, field.nNodesThreshold
#   xEgo, yEgo, psiEgo, vEgo = zCurr

#   # Polynomial fit for each lane
#   coeffsLanes = Inf*ones(4,nLanes)
#   idxClosestPointLanes = round(Int64,zeros(nLanes))
#   for n=1:nLanes
#     # nodes for each lane
#     nodesLane = nodes[:,wayPointers[1,n]:wayPointers[2,n]]
#     nNodesLane = size(nodesLane,2)

#     # distances from current position
#     dists = sqrt(sum((nodesLane-repmat([xEgo,yEgo],1,nNodesLane)).^2,1))

#     # distance and index of closest nodes
#     minDist, minIdx = findmin(dists)
#     idxClosestPointLanes[n] = wayPointers[1,n] + minIdx - 1
#     minIdx = min(minIdx, nNodesLane-1)
#     minIdx = max(minIdx, 2)

#     # direction of motion
#     dirEgo = [cos(psiEgo), sin(psiEgo)]
#     dirPos = nodesLane[:,minIdx+1] - nodesLane[:,minIdx]
#     direction = (dirEgo'*dirPos)[1] < 0 ? -1 : 1

#     # nodes for polynomial fit
#     idxStart = max(1, minIdx - (direction == 1 ? nNodesPolyBack : nNodesPolyFront))
#     idxEnd = min(nNodesLane, minIdx + (direction == 1 ? nNodesPolyFront : nNodesPolyBack))
#     nodesNear = nodesLane[:,idxStart:idxEnd]
#     nNodesNear = size(nodesNear,2)

#     if nNodesNear >= nNodesThreshold
#       # transform points to ego frame
#       R = [cos(psiEgo) -sin(psiEgo); sin(psiEgo) cos(psiEgo)]
#       nodesLocal = R'*(nodesNear - repmat([xEgo,yEgo],1,nNodesNear))

#       # compute least squares fit
#       px = nodesLocal[1,:]'
#       py = nodesLocal[2,:]'
#       H = [ones(nNodesNear) px px.^2 px.^3]
#       coeffsLanes[:,n] = -(H'*H)\H'*py
#     end
#   end

#   # closest lane index
#   minEy, laneIdx = findmin(abs(coeffsLanes[1,:]))
#   coeffsCurrLane = coeffsLanes[:,laneIdx]
#   mapSpeedRef = speed[idxClosestPointLanes[laneIdx]]

#   return coeffsCurrLane, laneIdx, mapSpeedRef
# end

# ### Generate state and input reference trajectory for MPC problem
# function generateReference(zCurr::Array{Float64,1},field::TrackMap,
#                            tuning::Tuning,mpc::MpcProblem)
#   # parameters
#   dt, N = tuning.dt, tuning.N
#   xEgo, yEgo, psiEgo, vEgo = zCurr
#   nz, nu = mpc.nz, mpc.nu
#   xLocal = field.xLocal

#   # localize vehicle in lane and get reference speed
#   coeffsCurrLane, laneIdx, vRef = localizeVehicle(zCurr,field)

#   # longitudinal distance along lane centerline as a function of X
#   yLocal = polyval(Poly(-coeffsCurrLane),xLocal)
#   sLocal = [0.0; cumsum(sqrt(diff(xLocal).^2 + diff(yLocal).^2))]

#   # desired longitudinal positions
#   sDesired = [0.0; cumsum(dt*vRef*ones(N))]

#   # reference points on lane centerline
#   itpX = interpolate((sLocal,), xLocal, Gridded(Linear()))
#   itpY = interpolate((sLocal,), yLocal, Gridded(Linear()))
#   xLaneLocal, yLaneLocal = itpX[sDesired], itpY[sDesired]
#   posRefLocal = hcat(xLaneLocal,yLaneLocal)'

#   # reference points in global frame
#   R = [cos(psiEgo) -sin(psiEgo); sin(psiEgo) cos(psiEgo)]
#   posRefGlobal = repmat([xEgo,yEgo],1,N+1) + R*posRefLocal

#   # heading reference
#   psiRef = psiEgo*ones(N+1,1) + [0.0; atan2(diff(xLaneLocal),diff(yLaneLocal))]

#   # state reference
#   ZRef = [posRefGlobal; psiRef'; vRef*ones(1,N+1)]
#   URef = zeros(nu,N)

#   # outputs
#   return URef, ZRef
# end

end # module
