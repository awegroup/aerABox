% This function generates an APAME geometry and the input file given
% certain parameters
%       Author  : Gabriel Buendia
%       Version : 1
% Inputs:
%       Slices_Wing -> Matrix [N_points_airfoil x 3 x N_slices] containing 
%                      all the wing cross sections [m]
%       h           -> Height between wings [m]
%       name        -> Name of the input file
%       airspeed    -> Freestream reference speed [m/s]
%       density     -> Freestream density [kg/m^3]
%       pressure    -> Freestream pressure [Pa]
%       mach        -> Fresstream Mach number
%       cases       -> Vector [n x 2] of cases to simulate in APAME 
%                      containing the AoA [deg] and AoS [deg] to simulate
%       wingspan    -> Wingspan [m]
%       MAC         -> Mean aerodynamic chord [m]
%       surf        -> Surface area [m]
%       origin      -> Origin of coordinates [m]
%       method      -> # singularity method:
%                      0-constant source/doublet 
%                      1-constant doublet
%       err         -> Error
%       colldepth   -> Collocation point depth [m] 
%       farfield    -> Far field coefficient
%       collcalc    -> Collocation point calculation:
%                      0-approximate
%                      1-accurate
%       velorder    -> Interpolation method/order for velocity calculations:
%                      0-nodal
%                      1-first
%                      2-second
%       results     -> Result requests:
%                      0-no
%                      1-yes
%       requests    -> Type of request: [vector 1x13] specify 1 or 0
%                      1-coefficients
%                      2-forces
%                      3-geometry
%                      4-velocity
%                      5-pressure
%                      6-center points
%                      7-dipole values 
%                      8-source values
%                      9-velocity components
%                      10-mesh characteristics
%                      11-static pressure
%                      12-dynamic pressure 
%                      13-manometer pressure
% Outputs:
%        Nopan -> Number of wing panels
function [NoPan,area] = BoxWingAPAME(Slices_Wing,h,name,airspeed,density,pressure,mach,cases,wingspan,MAC,surf,origin,method,err,colldepth,farfield,collcalc,velorder,results,requests)
addpath("sr_bwa\")

% Wing Panels
[Nx, Ny, Nz] = size(Slices_Wing);
for i = 1:Nz
    Slices_Wing(:,:,i) = Slices_Wing(linspace(Nx,1,Nx),:,i);
end
Slices_Wing(:,3,:) = Slices_Wing(:,3,:)-h/2;
aux_z = Nz;
for k = 1:Nz
    aux_z = aux_z + 1;
    for j = 1:Ny
        for i = 1:Nx
            if j == 2
                Slices_Wing(i,j,aux_z) = -Slices_Wing(i,j,Nz-k+1);
            else
                Slices_Wing(i,j,aux_z) = Slices_Wing(i,j,Nz-k+1);
            end
        end
    end
end
[Nx, Ny, Nz] = size(Slices_Wing);
aux_nodes = 0;
Wake_Points = 5;
Spacing_Wake = 10;
Nodes = (Nx-1)*Nz;
Nodes_List = zeros(Nodes,3);
Track_Nodes = zeros(Nx-1,Nz);
Track_Panels = zeros(Nx-1,Nz-2);
area = 0;
for k = 1:Nz
    for i = 1:Nx-1
        aux_nodes = aux_nodes+1;
        Nodes_List(aux_nodes,:) = Slices_Wing(i,:,k);
        Track_Nodes(i,k) = aux_nodes;
    end
end
aux_panels = 0;
for j = 1:Nz-2
    for i = 1:Nx-1
        aux_panels = aux_panels + 1;
        Track_Panels(i,j) = aux_panels;
    end
end
aux_pannodes = 0;
Panels = (Nz-2)*(Nx-1);
NoPan = Panels;
Panel_Nodes = zeros(Panels,4);
for j = 1:Nz-2
    for i = 1:Nx-1
        aux_pannodes = aux_pannodes + 1;
        if j < Nz/2
            if i == Nx-1
               Panel_Nodes(aux_pannodes,1) = Track_Nodes(1,j);
               Panel_Nodes(aux_pannodes,4) = Track_Nodes(1,j+1);  
               Panel_Nodes(aux_pannodes,2) = Track_Nodes(i,j);
               Panel_Nodes(aux_pannodes,3) = Track_Nodes(i,j+1);
            elseif i > Nx/2
               Panel_Nodes(aux_pannodes,1) = Track_Nodes(i+1,j);
               Panel_Nodes(aux_pannodes,4) = Track_Nodes(i+1,j+1);
               Panel_Nodes(aux_pannodes,2) = Track_Nodes(i,j);
               Panel_Nodes(aux_pannodes,3) = Track_Nodes(i,j+1);
            else
               Panel_Nodes(aux_pannodes,4) = Track_Nodes(i+1,j);
               Panel_Nodes(aux_pannodes,3) = Track_Nodes(i+1,j+1);
               Panel_Nodes(aux_pannodes,1) = Track_Nodes(i,j);
               Panel_Nodes(aux_pannodes,2) = Track_Nodes(i,j+1); 
            end
        elseif j == Nz-2
            if i == Nx-1
                Panel_Nodes(aux_pannodes,1) = Track_Nodes(1,j+1);
                Panel_Nodes(aux_pannodes,4) = Track_Nodes(1,end);
                Panel_Nodes(aux_pannodes,2) = Track_Nodes(i,j+1);
                Panel_Nodes(aux_pannodes,3) = Track_Nodes(i,end);
            elseif i > Nx/2
                Panel_Nodes(aux_pannodes,1) = Track_Nodes(i+1,j+1);
                Panel_Nodes(aux_pannodes,4) = Track_Nodes(i+1,end);    
                Panel_Nodes(aux_pannodes,2) = Track_Nodes(i,j+1);
                Panel_Nodes(aux_pannodes,3) = Track_Nodes(i,end);
            else
                Panel_Nodes(aux_pannodes,2) = Track_Nodes(i+1,j+1);
                Panel_Nodes(aux_pannodes,1) = Track_Nodes(i+1,end);    
                Panel_Nodes(aux_pannodes,3) = Track_Nodes(i,j+1);
                Panel_Nodes(aux_pannodes,4) = Track_Nodes(i,end);
            end
        else
            if i == Nx-1
               Panel_Nodes(aux_pannodes,1) = Track_Nodes(1,j+1);
               Panel_Nodes(aux_pannodes,4) = Track_Nodes(1,j+2);  
               Panel_Nodes(aux_pannodes,2) = Track_Nodes(i,j+1);
               Panel_Nodes(aux_pannodes,3) = Track_Nodes(i,j+2);
            elseif i > Nx/2
               Panel_Nodes(aux_pannodes,1) = Track_Nodes(i+1,j+1);
               Panel_Nodes(aux_pannodes,4) = Track_Nodes(i+1,j+2);
               Panel_Nodes(aux_pannodes,2) = Track_Nodes(i,j+1);
               Panel_Nodes(aux_pannodes,3) = Track_Nodes(i,j+2);
            else
               Panel_Nodes(aux_pannodes,2) = Track_Nodes(i+1,j+1);
               Panel_Nodes(aux_pannodes,1) = Track_Nodes(i+1,j+2);
               Panel_Nodes(aux_pannodes,3) = Track_Nodes(i,j+1);
               Panel_Nodes(aux_pannodes,4) = Track_Nodes(i,j+2); 
            end
        end
        area_aux = quadrilateralArea(Nodes_List(Panel_Nodes(aux_pannodes,:),:));
        area = area + area_aux;
    end
end
aux_panconnect = 0;
Panel_Connectivity = zeros(Panels,4);
for j = 1:Nz-2
    for i = 1:Nx-1
        aux_panconnect = aux_panconnect + 1;
        if i == 1 
            Panel_Connectivity(aux_panconnect,4) = 0;
            Panel_Connectivity(aux_panconnect,2) = Track_Panels(i+1,j);
            if j == 1
                Panel_Connectivity(aux_panconnect,1) = Track_Panels(i,j+1);
                Panel_Connectivity(aux_panconnect,3) = Track_Panels(i,end);
            elseif j == Nz-2
                Panel_Connectivity(aux_panconnect,1) = Track_Panels(i,1);
                Panel_Connectivity(aux_panconnect,3) = Track_Panels(i,j-1);
            else
                Panel_Connectivity(aux_panconnect,1) = Track_Panels(i,j+1);
                Panel_Connectivity(aux_panconnect,3) = Track_Panels(i,j-1);
            end
        elseif i == Nx-1
            Panel_Connectivity(aux_panconnect,1) = Track_Panels(i-1,j);
            Panel_Connectivity(aux_panconnect,4) = 0;
            if j == 1
                Panel_Connectivity(aux_panconnect,2) = Track_Panels(i,j+1);
                Panel_Connectivity(aux_panconnect,3) = Track_Panels(i,end);
            elseif j == Nz-2
                Panel_Connectivity(aux_panconnect,2) = Track_Panels(i,1);
                Panel_Connectivity(aux_panconnect,3) = Track_Panels(i,j-1);
            else
                Panel_Connectivity(aux_panconnect,2) = Track_Panels(i,j+1);
                Panel_Connectivity(aux_panconnect,3) = Track_Panels(i,j-1);
            end
        else
            Panel_Connectivity(aux_panconnect,1) = Track_Panels(i-1,j);
            Panel_Connectivity(aux_panconnect,3) = Track_Panels(i+1,j);
            if j == 1
                Panel_Connectivity(aux_panconnect,2) = Track_Panels(i,j+1);
                Panel_Connectivity(aux_panconnect,4) = Track_Panels(i,end);
            elseif j == Nz-2
                Panel_Connectivity(aux_panconnect,2) = Track_Panels(i,1);
                Panel_Connectivity(aux_panconnect,4) = Track_Panels(i,j-1);
            else
                Panel_Connectivity(aux_panconnect,2) = Track_Panels(i,j+1);
                Panel_Connectivity(aux_panconnect,4) = Track_Panels(i,j-1);
            end
        end
    end
end
Panel_List = ones(Panels,9);

% Wake panels
Panels = Panels + Wake_Points*(Nz-2);
Nodes  = Nodes + Wake_Points*Nz;

for k = 1:Nz
    for i = 1:Wake_Points
        aux_nodes = aux_nodes+1;
        v_aux1 = Slices_Wing(1,:,k)-Slices_Wing(2,:,k); v_aux1 = v_aux1/norm(v_aux1);
        v_aux2 = Slices_Wing(end,:,k)-Slices_Wing(end-1,:,k); v_aux2 = v_aux2/norm(v_aux2);
        v_aux_avg = (v_aux1+v_aux2)/2; v_aux_avg = v_aux_avg/norm(v_aux_avg);
        aux_wake = Spacing_Wake*i*v_aux_avg;
        Nodes_List(aux_nodes,:) = aux_wake+Slices_Wing(1,:,k);
        Track_Nodes(Nx-1+i,k) = aux_nodes;
    end
end

for j = 1:Nz-2
    for i = Nx-1+1:Nx-1+Wake_Points
        aux_pannodes = aux_pannodes + 1;
        if j < Nz/2
            if i == Nx-1+1
               Panel_Nodes(aux_pannodes,2) = Track_Nodes(1,j);
               Panel_Nodes(aux_pannodes,1) = Track_Nodes(1,j+1);  
               Panel_Nodes(aux_pannodes,3) = Track_Nodes(i,j);
               Panel_Nodes(aux_pannodes,4) = Track_Nodes(i,j+1);
            else
               Panel_Nodes(aux_pannodes,3) = Track_Nodes(i,j);
               Panel_Nodes(aux_pannodes,4) = Track_Nodes(i,j+1);
               Panel_Nodes(aux_pannodes,2) = Track_Nodes(i-1,j);
               Panel_Nodes(aux_pannodes,1) = Track_Nodes(i-1,j+1); 
            end
        elseif j == Nz-2
            if i == Nx-1+1
                Panel_Nodes(aux_pannodes,4) = Track_Nodes(1,j+1);
                Panel_Nodes(aux_pannodes,3) = Track_Nodes(1,end);
                Panel_Nodes(aux_pannodes,1) = Track_Nodes(i,j+1);
                Panel_Nodes(aux_pannodes,2) = Track_Nodes(i,end);
            else
                Panel_Nodes(aux_pannodes,1) = Track_Nodes(i,j+1);
                Panel_Nodes(aux_pannodes,2) = Track_Nodes(i,end);    
                Panel_Nodes(aux_pannodes,4) = Track_Nodes(i-1,j+1);
                Panel_Nodes(aux_pannodes,3) = Track_Nodes(i-1,end);
            end
        else
            if i == Nx-1+1
               Panel_Nodes(aux_pannodes,4) = Track_Nodes(1,j+1);
               Panel_Nodes(aux_pannodes,3) = Track_Nodes(1,j+2);  
               Panel_Nodes(aux_pannodes,1) = Track_Nodes(i,j+1);
               Panel_Nodes(aux_pannodes,2) = Track_Nodes(i,j+2);
            else
               Panel_Nodes(aux_pannodes,1) = Track_Nodes(i,j+1);
               Panel_Nodes(aux_pannodes,2) = Track_Nodes(i,j+2);
               Panel_Nodes(aux_pannodes,4) = Track_Nodes(i-1,j+1);
               Panel_Nodes(aux_pannodes,3) = Track_Nodes(i-1,j+2); 
            end
        end
    end
end
for j = 1:Nz-2
    for i = Nx-1+1:Nx-1+Wake_Points
        aux_panconnect = aux_panconnect + 1;
        Panel_Connectivity(aux_panconnect,1) = Track_Panels(1,j);
        Panel_Connectivity(aux_panconnect,2) = Track_Panels(Nx-1,j);
        Panel_Connectivity(aux_panconnect,3) = NaN;
        Panel_Connectivity(aux_panconnect,4) = NaN;
    end
end

for i = 1:Panels
    Panel_List(i,2:5) = Panel_Nodes(i,:);
    Panel_List(i,6:9) = Panel_Connectivity(i,:);
    if Panel_List(i,1) == 0
        Panel_List(i,1) = 10;
    end
end
commandAPAME = inputFileAPAME(airspeed,density,pressure,mach,cases,wingspan,MAC,surf,origin,method,err,colldepth,farfield,collcalc,velorder,results,requests,Nodes,Nodes_List,Panels,Panel_List);

% Writing in file
fileName = fopen(name,'w');
fprintf(fileName,commandAPAME);
fclose(fileName);
end