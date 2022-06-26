clear; close all; clc;
% APAME path
pathAPAME = '"F:\Universidad\Master_in_Aerospace_Engineering\TFM\Apame_v140915\ApameSolver\bin\apame_win64.exe"';
start = -14.8;
finish = 9.2;
chord = 3.5425;
surf_lim = 150.4478;
span = 21.2347;
sweep_b =linspace(span,2*span,16);
vec = linspace(start,finish,3);
data_APAME = zeros(length(vec),4,length(sweep_b));
for kindex = 1:length(sweep_b)
close all
b = sweep_b(kindex);
h = 0.2*b;                                    
fore_area = 0.5;
taper = 1;
c_r1 = fore_area*surf_lim*2/(b*(1+taper));                                  
c_r2 = (1-fore_area)*surf_lim*2/(b*(1+taper));             
c_t1 = taper*c_r1;                              
c_t2 = taper*c_r2;                                
N_bu = 3;                                  
N_bl = 3;                                   
N_bw = 0;                                   
Finite = 0;                                 
N_el = 30;
N = 5;
af_1r = '23012';
af_2r = '23012';
af_1t = '23012';
af_2t = '23012';
af_w = '0030';

AoS = 0;
height = 0;
[T,~,p,~] = atmosisa(height);
[rho,mu] = AirProperties(T-273.85, p/100, [],'rho', 'mu');
gamma = 1.4;
R_c = 287.04;
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

s = 0.3*MAC;
R = 0.15*sqrt(h^2+s^2);                                  

airspeed = 60;
Re = rho*MAC*airspeed/mu;
mach = airspeed/sqrt(gamma*R_c*T);
origin = [0 0 0];
method = 0;
err = 1e-7;
colldepth = 1e-7;
farfield = 5;
collcalc = 1;
velorder = 2;
results = 1;
requests = zeros(1,13); requests(1) = 1;

[Slices_Wing, Points_Wing] = ParametrizationBWFiniteTE(h,b,R,c_r1,c_r2,c_t1,c_t2,s,N_bu,N_bl,N_bw,N_el,N,Finite,af_1r,af_2r,af_1t,af_2t,af_w);

counter = 0;
for indAoA = vec
    counter = counter + 1;
    AoA = indAoA;
    if indAoA < 0
        name = append('BWAPAME_High_n',num2str(abs(indAoA)),'.inp');
        APAME = append('BWAPAME_High_n',num2str(abs(indAoA)));
        string = append('BWAPAME_High_n',num2str(abs(indAoA)),'.res');
    else
        name = append('BWAPAME_High_',num2str(indAoA),'.inp');
        APAME = append('BWAPAME_High_',num2str(indAoA));
        string = append('BWAPAME_High_',num2str(indAoA),'.res');
    end
    cases = [AoA; AoS];
    [NoPan,area] = BoxWingAPAME(Slices_Wing,h,name,airspeed,density,pressure,mach,cases,wingspan,MAC,surf,origin,method,err,colldepth,farfield,collcalc,velorder,results,requests);
    cmdname = 'command.txt';
    cmdName = fopen(cmdname,'w');
    fprintf(cmdName,APAME);
    fclose(cmdName);
    system(append(pathAPAME,' < command.txt'));
    fileID = fopen(string,'r');
    string_full = fscanf(fileID,'%s');
    fclose(fileID);
    index_CD = strfind(string_full,'CDRAG');
    index_CD = index_CD+5;
    value_CD = string_full(index_CD:index_CD+13);
    value_CD = str2double(value_CD);
    index_CL = strfind(string_full,'CLIFT');
    index_CL = index_CL+5;
    value_CL = string_full(index_CL:index_CL+13);
    value_CL = str2double(value_CL);
    data_APAME(counter,1,kindex) = AoA; data_APAME(counter,2,kindex) = value_CL; data_APAME(counter,3,kindex) = value_CD; data_APAME(counter,4,kindex) = NoPan;
end
end
[C_L_alpha,C_L_0,valid] = liftCurve(data_APAME);
save('EquivalentWingMEGAWESSweepb.mat')

%% Plots
clear; close all; clc;
load('MEGAWESData.mat')
load('EquivalentWingMEGAWESSweepb.mat')
figure()
hold on
plot(alpha,cL,'ok','LineWidth',2)
cL_a = (cL(end)-cL(1))/((finish-start)*pi/180);
cL_0 = interp1(alpha,cL,0);
vector = {'-b','--b','ob','-r','--r','or','-g','--g','og','-m','--m','om','-c','--c','oc','-y','--y','oy'};
aux_c = 0;
for j = 1:kindex
    C_L_APAME = data_APAME(:,2,j);
    alpha_p = data_APAME(:,1,j);
    if valid(j) == 1
        aux_c = aux_c+1;
        plot(alpha_p,C_L_APAME,vector{j},'LineWidth',2)
        aux_sweep_b(aux_c) = sweep_b(j);
        [alpha_I,cL_I] = intersectionCurves(cL_a,cL_0,C_L_alpha(j),C_L_0(j));
        Int(aux_c) = liftDiffIntegral(alpha_I,cL,C_L_APAME,alpha_p*pi/180);
    end
    C_L_APAME = 0;
    alpha_p = 0;
end
xlabel('\alpha','FontSize',14)
ylabel('C_L','FontSize',14)
set(gca,'FontSize',14)
grid on
figure()
plot(aux_sweep_b,Int,'o','LineWidth',2)
xlabel('b [m]','FontSize',14)
ylabel('I','FontSize',14)
set(gca,'FontSize',14)
grid on
figure()
yyaxis left
plot(sweep_b,C_L_alpha(:),'o','LineWidth',2)
hold on 
plot(2*span,cL_a,'x','LineWidth',2)
ylabel('C_{L\alpha}','FontSize',14)
yyaxis right
plot(sweep_b,C_L_0(:),'o','LineWidth',2)
hold on 
plot(2*span,cL_0,'x','LineWidth',2)
xlabel('b [m]','FontSize',14)
ylabel('C_{L0}','FontSize',14)
set(gca,'FontSize',14)
grid on

function [C_L_alpha,C_L_0,valid] = liftCurve(data_APAME)
[a,b,c,d] = size(data_APAME);
for i = 1:c
    for j = 1:d
        C_L_APAME = data_APAME(:,2,i,j);
        alpha_APAME = data_APAME(:,1,i,j);
        [~,R1]=fit(alpha_APAME,C_L_APAME,'poly1');
        [~,R2]=fit(alpha_APAME,C_L_APAME,'poly2');
        C_L_alpha(i,j) = (C_L_APAME(end)-C_L_APAME(1))/((alpha_APAME(end)-alpha_APAME(1))*pi/180);
        C_L_0(i,j) = interp1(alpha_APAME,C_L_APAME,0);
        if R1.rsquare > 0.9995 && C_L_alpha(i,j) < 10 && C_L_alpha(i) >= 0
            valid(i,j) = 1;
        else
            valid(i,j) = 0;
        end
    end
end
end
function [x_I,y_I] = intersectionCurves(m1,n1,m2,n2)
    x_I = (n1-n2)/(m2-m1);
    y_I = m1*x_I+n1;
end
function I = liftDiffIntegral(alpha_I,cL,cL_A,alpha)

if alpha_I >= alpha(1) && alpha_I <= alpha(end) 
    I = (alpha_I-alpha(1))*(cL_A(1)-cL(6))/2 + (alpha(end)-alpha_I)*(cL_A(end)-cL(30))/2;
else
    I = ((cL_A(1)-cL(6))+(cL_A(end)-cL(30)))*(alpha(end)-alpha(1))/2;
end
end