% This function swaps the elements of a matrix
%       Author  : Gabriel Buendia
%       Version : 1
% Inputs:
%       A -> Matrix to swap
% Outputs:
%       A -> Swapped matrix
function A = swapmatrix(A)
[a,~] = size(A);
b = mod(a,2);
if b == 0
    a = a/2;
else
    a = (a-1)/2;
end
for i = 1:a
    aux = A(i,:);
    A(i,:) = A(end-(i-1),:);
    A(end-(i-1),:) = aux;
end