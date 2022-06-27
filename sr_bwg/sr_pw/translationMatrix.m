% This function aims to generate a translation matrix for a rigid body in
% cartesian coordinates
%       Author  : Gabriel Buendia
%       Version : 1
% Inputs:
%       tx: Translation vector component in x axis
%       ty: Translation vector component in y axis
%       tz: Translation vector component in z axis
% Outputs:
%       T : Translation matrix -> 3x3 Matrix 
function T = translationMatrix(tx,ty,tz)
T = [1 0 0 tx;
     0 1 0 ty;
     0 0 1 tz;
     0 0 0 1];
end