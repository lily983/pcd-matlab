function quat = getRandomQuaternion
% getRandomQuaternion Generate random quaternion for rotation
quat=rand(1,4);
quat = quat/norm(quat);
end