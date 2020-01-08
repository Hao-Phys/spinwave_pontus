function [qx, qy] = q12toqxy(q1, q2)
% this function is a transformation of change of basis
% from the basis of S-S lattice to the global reference frame

qx = 2.0*pi*q1;
qy = 4.0*pi/sqrt(3)*(q2-0.5*q1);
end