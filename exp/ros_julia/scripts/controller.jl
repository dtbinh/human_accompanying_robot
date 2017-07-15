#!/usr/bin/env julia

# controller for the robot
# originally from MPC
# Chang Liu 6/7/16
# modified for ROS use
# Chang Liu 6/30/17

using RobotOS
@rosimport std_msgs.msg: Float64MultiArray, Float64Msg
rostypegen()
using std_msgs.msg

using CarMpcUtils
using CarMpc

# ################
# ##### Setup #####
# ################

### Simulation parameters
simLength = 5
# initialize vehicle
robot = Robot([1.0, 1.0, 0.0, 0.5], model=unicycleDiscrete)

### tuning parameters
tuning = Tuning(dt = 0.5, N = 10, safe_dis = 1.0, safe_margin = 1.0, cmft_dis = 2.4)

### map
field = Field("maps/RFS_2Lanes_Speed_0128.mat")

### obstacles
obs = Array{Obstacle}(2,1)
obs[1] = Obstacle([3.0,5.0,0.0,0.0],"circle",0.0)
obs[2] = Obstacle([2.0,2.0,0.0,0.0],"circle",0.0)

### Initialize MPC problem and solve dummy problem
mpc = initializeMPCProblem(robot, field, tuning)
nz, nu, N = mpc.nz, mpc.nu, tuning.N

### Subscriber callbacks
callback(msg::Float64MultiArray) = begin
    init_state = msg.data;
end

sub_state_callback(sub_msg::Float64MultiArray, z0::Array{Float64,1}) = begin
  z0[1:nz] = sub_msg.data
end

################
##### Main #####
################

### MPC model parameters updated by subscribers
z0 = robot.z
u0 = zeros(nu)
z_h = zeros(4,N+1) ## update by ROS

USim = zeros(nu,simLength)
ZSim = [z0 zeros(nz,simLength)]

### ROS subroutines
init_node("mpc_julia")
pub_steer = Publisher{Float64MultiArray}("mpc_out_steer", queue_size=1)
pub_acc = Publisher{Float64MultiArray}("mpc_out_acc", queue_size=1)
sub = Subscriber{Float64MultiArray}("mpc_solver_feed", sub_state_callback, (z0,z_h), queue_size=1)

### Main loop
t = 1;
while !is_shutdown()
  ### Update and solve MPC problem
  # ZOpt, UOpt, solveTime = updateSolveMpcProblem(mpc, z0, u0)
  ZOpt, UOpt, solveTime = updateSolveMpcProblem(mpc, z0, z_h)

  ### Publish outputs
  steer_msg = Float64Msg()
  steer_msg.data = UOpt[1,1]
  publish(pub_steer, steer_msg)

  acc_msg = Float64Msg()
  acc_msg.data = UOpt[2,1]
  publish(pub_acc, acc_msg)

  ### Update current input
  u0[1:nu] = UOpt[:,1]

  ### Variables for logging
  ZSim[:,t+1] = robot.z
  USim[:,t] = u0
  print(robot.z)
  t = t+1;
end
