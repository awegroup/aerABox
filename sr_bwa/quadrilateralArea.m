% This function computes the area of a quadrilateral element given 4
% points
%       Author  : Gabriel Buendia
%       Version : 1
% Inputs:
%       points -> Vector [4x3] containing the coordinates [m] of the four
%                 points defining the quadrilateral element
% Outputs:
%       area -> Area of the quadrilateral element [m^2]
function area = quadrilateralArea(points)
vector = zeros(4,3);
vector(1,:) = points(2,:)-points(1,:);
vector(2,:) = points(4,:)-points(1,:);
vector(3,:) = points(4,:)-points(3,:);
vector(4,:) = points(2,:)-points(3,:);
cross_vec(1,:) = cross(vector(1,:),vector(2,:));
cross_vec(2,:) = cross(vector(3,:),vector(4,:));
area_1 = 0.5*norm(cross_vec(1,:));
area_2 = 0.5*norm(cross_vec(2,:));
area = area_1 + area_2;
end