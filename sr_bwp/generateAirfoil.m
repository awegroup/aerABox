% This function generates a 4 or 5 series NACA airfoil
%       Author  : Gabriel Buendia
%       Version : 1
% Inputs:
%       designation -> Airfoil designation
%       digits      -> 4 or 5
%       panels      -> Number of points discretizing the airfoil
%       TE          -> Trailing edge
% Outputs:
%       airfoil     -> vector [panels x 2] containing the x and z
%                      coordinates of the airfoil
function airfoil = generateAirfoil(designation,digits,panels,TE)
iaf.designation= designation;
iaf.n= panels;
iaf.HalfCosineSpacing=1;
iaf.wantFile=0;
iaf.datFilePath='./';
iaf.is_finiteTE=TE;
switch digits
    case 4
        af = naca4gen(iaf);
    case 5
        af = naca5gen(iaf);
end
airfoil = zeros(iaf.n*2+1,2);
airfoil(:,1) = af.x;
airfoil(:,2) = af.z;
end