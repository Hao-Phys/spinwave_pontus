function [q1, q2] = qxytoq12(qx, qy)
% this function is a transformation of change of basis
% from the global frame to the basis of S-S lattice
q1 = qx/(2.0*pi);
q2 = sqrt(3)/(4.0*pi)*qy + qx/(4.0*pi);
end