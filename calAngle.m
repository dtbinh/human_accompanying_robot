% 7/11/14
% The function calculates the angle (in radian) of a vector

function angle = calAngle (vec)
len = norm(vec);
x = vec(1);
y = vec(2);
if x >= 0
    angle = asin(y/len);
else
    angle = pi - asin(y/len);
end
end