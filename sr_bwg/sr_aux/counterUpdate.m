% This function aim to reset the counter for each airfoil section
%       Author  : Gabriel Buendia
%       Version : 1
% Inputs:
%       N       -> Number of points in a section
%       counter -> Actual value of the counter at the end of the section
% Output:
%       counter -> Updated counter to the next decimal place
function counter = counterUpdate(N,counter)
digits = ceil(log10(N));
counter = counter/10^digits;
counter = ceil(counter);
counter = counter*10^digits;
end