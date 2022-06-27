% This function aims to generate a rotation matrix for a rigid body in
% cartesian coordinates
%       Author  : Gabriel Buendia
%       Version : 1
% Inputs:
%       rx: Rotation angle in x axis [deg]
%       ry: Rotation angle in y axis [deg]
%       rz: Rotation angle in z axis [deg]
% Outputs:
%       R : Rotation matrix -> 3x3 Matrix 
function R = rotationMatrix(rx,ry,rz)
Rx = [1 0 0 0;
      0 cosd(rx) sind(rx) 0;
      0 -sind(rx) cosd(rx) 0;
      0 0 0 1];
Ry = [cosd(ry) 0 -sind(ry) 0;
      0 1 0 0;
      sind(ry) 0 cosd(ry) 0;
      0 0 0 1];
Rz = [cosd(rz) -sind(rz) 0 0;
      sind(rz) cosd(rz) 0 0;
      0 0 1 0;
      0 0 0 1];
R = Rx*Ry*Rz;
end