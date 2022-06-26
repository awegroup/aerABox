% This function computes the distance between nodes in the BL region of a curved 
% 2D mesh depending on the layer, angle of the curve, growth rate and total length
%       Author  : Gabriel Buendia
%       Version : 1
% Inputs:
%       angle  -> Angle separating the curved region [º]
%       l_init -> Initial boundary layer thickness [m]
%       GR     -> Growth rate
%       L      -> Radius of the curve = Length of the segment [m]
% Outputs:
%       d -> Vector [1xN_elements_curve] containing the evolution of the 
%            element sizes depending on the layer of the BL
function d = lineDivision(angle,l_ini,GR,L)
i = 0;
l = 0;
while l <= L
    i = i+1;
    l = l_ini*GR^i;
    x1 = L-l;
    y1 = 0;
    x2 = (L-l)*cosd(angle);
    y2 = (L-l)*sind(angle);
    d(i) = sqrt((x2-x1)^2+(y2-y1)^2);
end
end