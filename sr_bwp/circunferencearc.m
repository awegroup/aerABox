% Function that computes N points of the arc of the circunference 
% given 2 orthogonal points and the center (Note that it only works 
% if both points are orthogonal). Note that the arc ends are excluded.
%       Author  : Gabriel Buendia
%       Version : 1
% Inputs: 
%       rc -> Center of the circunference
%       r1 -> Point 1
%       r2 -> Point 2
%       N  -> Number of points
% Output:
%       CP -> Points on the circunference
function CP = circunferencearc(rc,r1,r2,N)
theta = linspace(0,pi/2,N);
theta = theta(2:end-1);
x = @(t) rc(1) + cos(t).*(r1(1)-rc(1)) + sin(t).*(r2(1)-rc(1));
y = @(t) rc(2) + cos(t).*(r1(2)-rc(2)) + sin(t).*(r2(2)-rc(2));
z = @(t) rc(3) + cos(t).*(r1(3)-rc(3)) + sin(t).*(r2(3)-rc(3));
CP(:,1) = x(theta);
CP(:,2) = y(theta);
CP(:,3) = z(theta);
end