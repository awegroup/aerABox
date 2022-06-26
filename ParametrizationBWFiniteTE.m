% This function generates the geometry of a boxwing given certain
% parameters and stores it in two different ways: a matrix containing all
% the cross sections defining the wing and a list of all points present
% in the wing
%       Author  : Gabriel Buendia
%       Version : 1
% Inputs:
%       h        -> Height between wings [m]
%       b        -> Wingspan [m]
%       R        -> Cant radius [m]
%       c_r1     -> Chord at the root of the lower wing [m]
%       c_t1     -> Chord at the tip of the lower wing [m]
%       c_r2     -> Chord at the root of the upper wing [m]
%       c_t2     -> Chord at the tip of the upper wing [m]
%       s        -> Wing stagger [m]
%       N_bu     -> Number of spanwise sections upper wing
%       N_bl     -> Number of spanwise sections lower wing
%       N_bw     -> Number of spanwise sections winglet
%       N_el     -> Number of elements to discretize airfoils
%       N        -> Number of sections to discretize the cant region
%       Finite   -> Finite trailing edge: 1-yes | 0-no
%       af_1r    -> Airfoil at the root of the lower wing 
%       af_1t    -> Airfoil at the tip of the lower wing 
%       af_2r    -> Airfoil at the root of the upper wing 
%       af_2t    -> Airfoil at the tip of the upper wing 
%       af_w     -> Airfoil at the winglet
% Outputs:
%       Slices_Wing -> Matrix [N_points_airfoil x 3 x N_slices] containing 
%                        all the wing cross sections [m]
%       Points_Wing -> Matrix [N_points_airfoil*N_slices x 3] containing
%                        all the points of the wing [m]
function [Slices_Wing, Points_Wing] = ParametrizationBWFiniteTE(h,b,R,c_r1,c_r2,c_t1,c_t2,s,N_bu,N_bl,N_bw,N_el,N,Finite,af_1r,af_2r,af_1t,af_2t,af_w)
addpath("sr_bwp\")

% Define the airfoils
% load("MEGAWESAirfoil.mat") % You can choose to load your own airfoil geometry, then you should specify it manually at each section
airfoil_1root = generateAirfoil(af_1r,length(af_1r),N_el,Finite);   % Dimensionless airfoil profile 1 root
% airfoil_1root = AirfoilRevEHC; % Replace your own airfoil variable name
airfoil_1root = swapmatrix(airfoil_1root);
airfoil_1tip = generateAirfoil(af_1t,length(af_1t),N_el,Finite);    % Dimensionless airfoil profile 1 tip
% airfoil_1tip = AirfoilRevEHC; % Replace your own airfoil variable name
airfoil_1tip = swapmatrix(airfoil_1tip);
airfoil_2root = generateAirfoil(af_2r,length(af_2r),N_el,Finite);  % Dimensionless airfoil profile 2 root
% airfoil_2root = AirfoilRevEHC; % Replace your own airfoil variable name
airfoil_2tip = generateAirfoil(af_2t,length(af_2t),N_el,Finite);    % Dimensionless airfoil profile 2 tip
% airfoil_2tip = AirfoilRevEHC; % Replace your own airfoil variable name
airfoil_winglet = generateAirfoil(af_w,length(af_w),N_el,Finite); % Dimensionless airfoil profile winglet
% airfoil_winglet = generateAirfoil(af_w,length(af_w),N_el,Finite); % Replace your own airfoil variable name

% Close open profiles (finite TE)
if Finite == 1
    prov_size = length(airfoil_winglet);
    airfoil_1root(prov_size+1,:) = airfoil_1root(1,:); 
    airfoil_1tip(prov_size+1,:) = airfoil_1tip(1,:); 
    airfoil_2root(prov_size+1,:) = airfoil_2root(1,:); 
    airfoil_2tip(prov_size+1,:) = airfoil_2tip(1,:); 
    airfoil_winglet(prov_size+1,:) = airfoil_winglet(1,:);
else
    airfoil_1root(1,2) = 0; airfoil_1root(end,2) = 0; 
    airfoil_1tip(1,2) = 0; airfoil_1tip(end,2) = 0; 
    airfoil_2root(1,2) = 0; airfoil_2root(end,2) = 0; 
    airfoil_2tip(1,2) = 0; airfoil_2tip(end,2) = 0; 
    airfoil_winglet(1,2) = 0; airfoil_winglet(end,2) = 0;
end

% Calculations
airfoil_1r = c_r1*airfoil_1root;            % Dimensional airfoil 1 at the root
airfoil_1t = c_t1*airfoil_1tip;             % Dimensional airfoil 1 at the tip
airfoil_2r = c_r2*airfoil_2root;            % Dimensional airfoil 2 at the root
airfoil_2t = c_t2*airfoil_2tip;             % Dimensional airfoil 2 at the tip
airfoil_1w = c_t1*airfoil_winglet;          % Dimensional airfoil winglet
airfoil_2w = c_t2*airfoil_winglet;          % Dimensional airfoil winglet
vec_h = [s 0 h]; vec_h = vec_h/norm(vec_h); % Versor connecting airfoils 1 and 2
vec_hn = [vec_h(3) 0 -vec_h(1)];            % Versor normal to vec_h
vec_i = [1 0 0];                            % Versor from x axis
vec_j = [0 1 0];                            % Versor from y axis
vec_k = [0 0 1];                            % Versor from z axis
p_O = [0 0 0];                              % Origin
p_A = [0 b/2-R 0];                          % End of lower wing
p_B = p_A + R*vec_h;                        % Center of circunference 1
p_C = p_B + R*vec_j;                        % Lower part of winglet
p_D = p_C + (sqrt(h^2+s^2)-2*R)*vec_h;      % Upper part of winglet
p_E = p_D - R*vec_j;                        % Center of circumference 2
p_F = p_E + R*vec_h;                        % End of upper wing
p_G = p_F - (b/2-R)*vec_j;                  % Root of upper wing

airfoil_O = zeros(length(airfoil_1r),3);
airfoil_O(:,1) = airfoil_1r(:,1); airfoil_O(:,3) = airfoil_1r(:,2);

airfoil_A = zeros(length(airfoil_1t),3)+p_A;
airfoil_A(:,1) = airfoil_1t(:,1)+p_A(1); airfoil_A(:,3) = airfoil_1t(:,2)+p_A(3);

airfoil_C = zeros(length(airfoil_1t),3)+p_C;
airfoil_C(:,1) = airfoil_1w(:,1)+p_C(1); airfoil_C(:,2) = airfoil_1w(:,2)+p_C(2);

airfoil_D = zeros(length(airfoil_2t),3)+p_D;
airfoil_D(:,1) = airfoil_2w(:,1)+p_D(1); airfoil_D(:,2) = airfoil_2w(:,2)+p_C(2);

airfoil_F = zeros(length(airfoil_2t),3)+p_F;
airfoil_F(:,1) = airfoil_2t(:,1)+p_F(1); airfoil_F(:,3) = airfoil_2t(:,2)+p_F(3);

airfoil_G = zeros(length(airfoil_2t),3)+p_G;
airfoil_G(:,1) = airfoil_2r(:,1)+p_G(1); airfoil_G(:,3) = airfoil_2r(:,2)+p_G(3);

center_temp = airfoil_A(1,:)+R*vec_h;
for i = 1:length(airfoil_1t)
    opts = optimoptions('fsolve', 'TolFun', 1E-10, 'TolX', 1E-10);
    center = fsolve(@(r0) functioncenter(r0,airfoil_A(i,:),airfoil_C(i,:),vec_hn),center_temp,opts);
    center_temp = center;
    CP1(:,:,i) = circunferencearc(center,airfoil_A(i,:),airfoil_C(i,:),N);
end
CP1(:,:,end) = CP1(:,:,1);
center_temp = airfoil_F(1,:)-R*vec_h;
for i = 1:length(airfoil_2t)
    opts = optimoptions('fsolve', 'TolFun', 1E-10, 'TolX', 1E-10);
    center = fsolve(@(r0) functioncenter(r0,airfoil_D(i,:),airfoil_F(i,:),vec_hn),center_temp,opts);
    center_temp = center;
    CP2(:,:,i) = circunferencearc(center,airfoil_D(i,:),airfoil_F(i,:),N);
end
CP2(:,:,end) = CP2(:,:,1);

% Join points and organize in lines
N_lines = length(airfoil_1r);
N_points_lines = 2*(N-2)+6+N_bu+N_bw+N_bl;
Points_Wing = zeros(N_lines*N_points_lines,3);
k = 0;
% Optional figure, suppress in optimization applications
% Uncomment lines 129, 130, 157, 193-199 for visualization
% figure()
% hold on
for i = 1:N_lines
        counter = i+(N_points_lines-1)*k; 
        Points_Wing(counter,:) = airfoil_O(i,:);
        aux_low = (airfoil_A(i,:)-airfoil_O(i,:))/(N_bl+1);
        for j  = 1:N_bl
            Points_Wing(counter+j,:) = airfoil_O(i,:)+j*aux_low;
        end
        Points_Wing(counter+1+N_bl,:) = airfoil_A(i,:);
        for j = 1:(N-2)
            Points_Wing(counter+1+N_bl+j,:) = CP1(j,:,i);
        end
        Points_Wing(counter+N+N_bl,:) = airfoil_C(i,:);
        aux_winglet = (airfoil_D(i,:)-airfoil_C(i,:))/(N_bw+1);
        for j  = 1:N_bw
            Points_Wing(counter+N+N_bl+j,:) = airfoil_C(i,:)+j*aux_winglet;
        end
        Points_Wing(counter+N+N_bl+N_bw+1,:) = airfoil_D(i,:);
        for j = 1:(N-2)
            Points_Wing(counter+j+N+N_bl+N_bw+1,:) = CP2(j,:,i);
        end
       Points_Wing(counter+2*N+N_bl+N_bw,:) = airfoil_F(i,:);
       aux_up = (airfoil_G(i,:)-airfoil_F(i,:))/(N_bu+1);
       for j  = 1:N_bu
           Points_Wing(counter+2*N+N_bl+N_bw+j,:) = airfoil_F(i,:)+j*aux_up;
       end
       Points_Wing(counter+2*N+N_bl+N_bw+N_bu+1,:) = airfoil_G(i,:);
%        plot3(Points_Wing(counter:counter+2*N+N_bl+N_bw+N_bu+1,1),Points_Wing(counter:counter+2*N+N_bl+N_bw+N_bu+1,2),Points_Wing(counter:counter+2*N+N_bl+N_bw+N_bu+1,3),'-k','LineWidth',0.5)
       k = k+1;
end

% Join points and organize in slices
N_slices = 2*(N-2)+6+N_bu+N_bl+N_bw;
N_points_slices = length(airfoil_A);
Slices_Wing = zeros(N_points_slices,3,N_slices);
k = 0;
Slices_Wing(:,:,1) = airfoil_O;
Increment_Low = (airfoil_A-airfoil_O)/(N_bl+1);
for i = 1:N_bl
    Slices_Wing(:,:,i+1) = airfoil_O+i*Increment_Low;
end
Slices_Wing(:,:,2+N_bl) = airfoil_A;
for i = 1:N-2
    aux = CP1(i,:,:);
    Slices_Wing(:,:,i+2+N_bl) = reshape(aux,[3 N_points_slices])';
end
Slices_Wing(:,:,N+1+N_bl) = airfoil_C;
Increment_Winglet = (airfoil_D-airfoil_C)/(N_bw+1);
for i = 1:N_bw
    Slices_Wing(:,:,i+N+1+N_bl) = airfoil_C+i*Increment_Winglet;
end
Slices_Wing(:,:,N+2+N_bl+N_bw) = airfoil_D;
for i = 1:N-2
    aux = CP2(i,:,:);
    Slices_Wing(:,:,N+i+2+N_bl+N_bw) =  reshape(aux,[3 N_points_slices])';
end
Slices_Wing(:,:,2*N+1+N_bl+N_bw) = airfoil_F;
Increment_Up = (airfoil_G-airfoil_F)/(N_bu+1);
for i = 1:N_bu
    Slices_Wing(:,:,i+2*N+1+N_bl+N_bw) = airfoil_F+i*Increment_Up;
end
Slices_Wing(:,:,2*N+2+N_bl+N_bw+N_bu) = airfoil_G;

% for i = 1:N_slices
%     plot3(Slices_Wing(:,1,i),Slices_Wing(:,2,i),Slices_Wing(:,3,i),'b','LineWidth',1.5)
% end
% grid on
% campos([-5.81744979365763,-7.348383400029182,6.328819187324961])
% drawnow
% axis equal