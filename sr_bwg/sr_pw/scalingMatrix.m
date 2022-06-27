% This function aims to generate a scaling matrix for a rigid body in
% cartesian coordinates.
%       Author  : Gabriel Buendia
%       Version : 1
% Inputs:
%       sx: Scaling factor in x axis 
%       sy: Scaling factor in y axis 
%       sz: Scaling factor in z axis 
% Outputs:
%       S : Scaling matrix -> 3x3 Matrix 
function S = scalingMatrix(sx,sy,sz)
S = [sx 0  0  0;
     0  sy 0  0;
     0  0  sz 0;
     0  0  0  1];
end