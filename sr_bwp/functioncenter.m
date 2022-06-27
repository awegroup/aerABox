% This function generates the system of equations to be solved for
% computing the center of rotation
%       Author  : Gabriel Buendia
%       Version : 1
% Inputs:
%       r0 -> Center of rotation
%       r1 -> Point 1 belonging to the arc
%       r2 -> Point 2 belonging to the arc
%       hn -> Vector normal to the circle plane
% Outputs:
%       F  -> System of five equations equated to 0
function F = functioncenter(r0,r1,r2,hn)
v1 = r1-r0;
v2 = r2-r0;
F(1) = dot(v1,v2);
F(2) = norm(v1)-norm(v2);
aux = cross(v1,v2)/norm(cross(v1,v2));
F(3) = aux(1)-hn(1);
F(4) = aux(2)-hn(2);
F(5) = aux(3)-hn(3);
end