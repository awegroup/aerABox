clear; close all; clc;
globalDirName = 'BoxWingPaper';
mkdir(globalDirName);
% Local paths
exportFolder = 'F:/Final_Code/'; % The bar sense is / for Pointwise
pathPointwise = '''F:\Programas\PointwiseV18.5R2\win64\bin''';
pathGlyphFile = 'F:\Final_Code\';
cases = 4;
for indexAlpha = 1:length(cases)
close all
factor = 0.0254;
h = 8*factor;                               % [m] Gap between the wings
b = 40*factor;                              % [m] Wingspan
R = 0.15*8*factor;                          % [m] Cant radius
c_r1 = 8*factor;                            % [m] Chord at the root of airfoil 1
c_r2 = 8*factor;                            % [m] Chord at the root of airfoil 2
c_t1 = 8*factor;                            % [m] Chord at the tip of airfoil 1
c_t2 = 8*factor;                            % [m]Chord at the tip of airfoil 2
s = 0;                                      % [m] Stagger
N_bu = 0;                                   % Number of spanwise sections upper wing
N_bl = 0;                                   % Number of spanwise sections lower wing
N_bw = 0;                                   % Number of spanwise sections winglet
Finite = 1;                                 % Finite TE
N_el = 40;                                  % Number of elements to discretize airfoils
N = 10;                                     % Number of sections to discretize the cant region
af_1r = '0012';
af_2r = '0012';
af_1t = '0012';
af_2t = '0012';
af_w = '0003';

AoA = cases(indexAlpha);
AoS = 0;
height = 0;
[T,~,p,~] = atmosisa(height);
[rho,mu] = AirProperties(T-273.85, p/100, [],'rho', 'mu');
gamma = 1.4;
R_c = 287.04;
Re = 510000;
density = rho;
pressure = p;
wingspan = b;
lambda_1 = c_t1/c_r1;
lambda_2 = c_t2/c_r2;
surf_1 = 2*(c_r1+c_t1)*b/4;
surf_2 = 2*(c_r2+c_t2)*b/4; 
surf = surf_1+surf_2;
MAC_1 = 2/3*c_r1*(1+lambda_1+lambda_1^2)/(1+lambda_1);
MAC_2 = 2/3*c_r2*(1+lambda_2+lambda_2^2)/(1+lambda_2);
MAC = (MAC_1*surf_1/surf+MAC_2*surf_2/surf);
airspeed = Re*mu/(rho*MAC);
mach = airspeed/sqrt(gamma*R_c*T);
Q = 0.5*rho*surf*airspeed^2;
nu = mu/rho;

[Slices_Wing, Points_Wing] = ParametrizationBWFiniteTE(h,b,R,c_r1,c_r2,c_t1,c_t2,s,N_bu,N_bl,N_bw,N_el,N,Finite,af_1r,af_2r,af_1t,af_2t,af_w);

% File
name = strcat('BoxWing.glf');

% Geometry Thresholds
tolBound = 1e-6;
tolThres = 1e-6;

% Boundary Layer
BLoptions.maxLayers = 34;
BLoptions.fullLayers = 34;
BLoptions.growthRate = 1.2;
BLoptions.push = 'true';
BLoptions.solverAttribute = 'AllAndConvertWallDoms';
BLoptions.firstLayerHeight = 8.69e-6;
cant_el_number = 1;
cant_el_size = pi*R/(2*(N-1)*cant_el_number);
cant_el_size_vec = lineDivision(90/((N-1)*cant_el_number),BLoptions.firstLayerHeight,BLoptions.growthRate,R);
el_size = 0.7*cant_el_size_vec(BLoptions.fullLayers);

% Surface Mesh Options
AR_cell = 1;
nodesTEz = 3;
advanced.IsoCellType = 'TriangleQuad';
advanced.Algorithm = 'AdvancingFrontOrtho';
advanced.EdgeMaximumLength = 'Boundary';
surfDom.wings = 'structured';
surfDom.cant = 'structured';
surfDom.winglet = 'structured';
distribution = 'Tanh';
refineVal = el_size/2;

% Farfield
FF.Type = 'Box';
FF.TypeCell = 'Unstructured';
FF.Dims = [40 15; 15 15; 15 15]*MAC;
symmetry = 'off';

% Isotropic Limits
ISOoptions.minLength = 'Boundary';
ISOoptions.maxLength = 1;
ISOoptions.maxGrowth = 1.2;

% Sources
Src1.Type = 'Box';
Src1.StartFace = 'x';
Src1.Center = [0.5*MAC 0 h/2];
Src1.Dims = [1.5*MAC 1.5*b 1.5*h];
Src1.Distr = 'Constant';
Src1.DistrVal = [0.004;0.3];
Src2.Type = 'Box';
Src2.StartFace = 'x';
rot_MAT = [cosd(AoA) 0 sind(AoA); 0 1 0;-sind(AoA) 0 cosd(AoA)];
Src2.Center = rot_MAT*Src1.Center';
Src2.Center = Src2.Center';
Src2.Center(1) = 5*MAC;
Src2.Dims = [10.5*MAC 1.5*b 2*h];
Src2.Distr = 'Parametric';
Src2.DistrVal = [0.0045 0.015;0.3 0.3];
sources = {Src1};

% Rotation
fixedPoint = [0 0 0];
rotAxis = 'y';
rotAngle = AoA;

% Export
string_AoA = num2str(AoA);
for i = 1:length(string_AoA)
    if isequal(string_AoA(i),'.') 
        string_AoA(i) = '_';
    end
    if isequal(string_AoA(i),'-')
        string_AoA(i) = 'n';
    end
end
mkdir(globalDirName,string_AoA);
save('Flight_Conditions.mat','Re','nu','p','rho','airspeed','surf','MAC')
exportPath = append(exportFolder,globalDirName,'/',string_AoA);
savePath = exportPath;
saveName = append('AoA_',string_AoA);
BoxWingGlyph(name,Slices_Wing,b,R,h,s,c_r1,c_t1,c_r2,c_t2,N,N_bl,N_bu,N_bw,el_size,tolBound,tolThres,AR_cell,nodesTEz,advanced,FF,symmetry,BLoptions,ISOoptions,sources,surfDom,distribution,refineVal,cant_el_number+1,fixedPoint,rotAxis,rotAngle,exportPath,savePath,saveName)
launchApplication = append('powershell -command "cd ',pathPointwise,'"; powershell -command .\Pointwise.exe ',pathGlyphFile,name);
system(launchApplication);
while ~exist(append(exportPath,'/',saveName,'.pw'))
    fprintf('...\n')
end
fprintf(append('Meshing completed and saved in ',exportPath,'/',saveName,'.pw\n'));
end