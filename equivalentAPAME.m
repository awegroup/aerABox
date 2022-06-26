clear; close all; clc;
% APAME path
pathAPAME = '"F:\Universidad\Master_in_Aerospace_Engineering\TFM\Apame_v140915\ApameSolver\bin\apame_win64.exe"';
% load('SweepBounded1.mat')
load('SweepBounded2.mat') % 1 Is the initial sweep (commented)
start = -14.8;
finish = 9.2;
AoA_vec = [start (start+finish)/2 finish];
surf_vec = zeros(length(combination),1);
area_vec = zeros(length(combination),1);
err_abs = zeros(length(combination),1);
data_APAME = zeros(length(AoA_vec),4,length(combination));
for i = 1:length(combination)
    counter = 0;
    b = combination(i,1); wingspan = b;
    MAC = combination(i,2); c_r1 = MAC; c_r2 = MAC; c_t1 = MAC; c_t2 = MAC;surf = 2*MAC*b;
    surf_vec(i) = surf;
    h = combination(i,3);
    s = combination(i,4);
    R = 0.15*sqrt(h^2+s^2);  
    N_bu = 3; N_bl = 3; N_bw = 0; N_el = 30; N = 5;                        
    Finite = 0;                               
    af_1r = '23012'; af_2r = '23012'; af_1t = '23012'; af_2t = '23012'; af_w = '0030';
    [Slices_Wing, Points_Wing] = ParametrizationBWFiniteTE(h,b,R,c_r1,c_r2,c_t1,c_t2,s,N_bu,N_bl,N_bw,N_el,N,Finite,af_1r,af_2r,af_1t,af_2t,af_w);
    AoS = 0; cases = [AoA; AoS];
    height = 0;
    [T,~,p,~] = atmosisa(height);
    [rho,mu] = AirProperties(T-273.85, p/100, [],'rho', 'mu');
    gamma = 1.4; R_c = 287.04; density = rho; pressure = p; 
    airspeed = 60; Re = rho*MAC*airspeed/mu; mach = airspeed/sqrt(gamma*R_c*T);
    origin = [0 0 0]; method = 0; err = 1e-7; colldepth = 1e-7; farfield = 5; collcalc = 1; velorder = 2; results = 1; requests = zeros(1,13); requests(1) = 1;    
    for indAoA = AoA_vec
        counter = counter +1;
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
        try
            index_CD = strfind(string_full,'CDRAG');
            index_CD = index_CD+5;
            value_CD = string_full(index_CD:index_CD+13);
            value_CD = str2double(value_CD);
            index_CL = strfind(string_full,'CLIFT');
            index_CL = index_CL+5;
            value_CL = string_full(index_CL:index_CL+13);
            value_CL = str2double(value_CL);
        catch
            value_CD = NaN;
            value_CL = NaN;
            AoA = NaN;
            NoPan = NaN;
        end
        data_APAME(counter,1,i) = AoA; data_APAME(counter,2,i) = value_CL; data_APAME(counter,3,i) = value_CD; data_APAME(counter,4,i) = NoPan;
    end
    area_vec(i) = area;
end
iaux = 0;
for i = 1:length(combination)
    if ~isnan(data_APAME(1,2,i))
        iaux = iaux + 1;
        aux_data_APAME(:,:,iaux) = data_APAME(:,:,i);
        aux_area_vec(iaux) = area_vec(i);
        aux_surf_vec(iaux) = surf_vec(i);
        aux_combination(iaux,:) = combination(i,:);
    end
end
save('EquivalentWingMEGAWESAirfoilFinal2.mat') % 1 Is the initial sweep (commented)

%% Plots
clear; close all; clc;
load('MEGAWESData.mat')
load('EquivalentWingMEGAWESAirfoilFinal2.mat') % 1 Is the initial sweep (commented)
data_APAME = aux_data_APAME;
area_vec = aux_area_vec;
surf_vec = aux_surf_vec;
filter_combination = aux_combination;
[C_L_alpha,C_L_0,valid] = liftCurve(data_APAME);

figure()
plot(alpha,cL,'ok','LineWidth',2)
hold on
cL_A = (cL(end)-cL(1))/((finish-start)*pi/180);
cL_0 = interp1(alpha*pi/180,cL,0);
vector = {'-b','--b','ob','-r','--r','or','-g','--g','og','-m','--m','om','-c','--c','oc','-y','--y','oy'};
aux_c = 0;
for i = 1:length(data_APAME(1,1,:))
        C_L_APAME = data_APAME(:,2,i);
        alpha_p = data_APAME(:,1,i);
        err_abs(i) = sum(abs(([cL(6) cL(18) cL(end)]-C_L_APAME')./[cL(6) cL(18) cL(end)]));
        h_b = filter_combination(i,3)/filter_combination(i,1);
        if valid(i) == 1 && h_b >= 0.1
            aux_c = aux_c+1;
            C_L_PANEL(aux_c,:) = C_L_APAME';
            plot(alpha_p,C_L_APAME,'-ob','LineWidth',2)
            aux_data_APAME2(:,:,aux_c) = data_APAME(:,:,i);
            aux_area_vec2(aux_c) = area_vec(i);
            aux_surf_vec2(aux_c) = surf_vec(i);
            aux_combination2(aux_c,:) = filter_combination(i,:);
            [alpha_I,cL_I] = intersectionCurves(cL_A,cL_0,C_L_alpha(i),C_L_0(i));
            Int(aux_c) = liftDiffIntegral(alpha_I,cL,C_L_APAME,alpha_p*pi/180);
        else
            C_L_alpha(i) = NaN;
            C_L_0(i) = NaN;
            err_abs(i) = NaN;
        end
        C_L_APAME = 0;
        alpha_p = 0;
end
data_APAME = aux_data_APAME2;
[C_L_alpha,C_L_0,valid] = liftCurve(data_APAME);
area_vec = aux_area_vec2;
surf_vec = aux_surf_vec2;
filter_combination = aux_combination2;
xlabel('\alpha','FontSize',14)
ylabel('C_L','FontSize',14)
set(gca,'FontSize',14)

n = 10;
[value_Max_1,index_Max_1] = nMaxValues(Int,length(filter_combination(:,1)));
count = 0;
diff_ini = 1;
diff_fin = 1;
while diff_ini > 0 && diff_fin > 0
    count = count + 1;
    diff_ini = data_APAME(1,2,index_Max_1(count))-cL(6);
    diff_fin = data_APAME(3,2,index_Max_1(count))-cL(30);
end
[value_Min_1,index_Min_1] = nMaxValues(1./area_vec(index_Max_1(1:count)),n);
 value_Min_1 = 1./value_Min_1;
figure(2)
plot(alpha,cL,'ok','LineWidth',2)
hold on
for i = 1:n
    figure(2)
    plot(AoA_vec,data_APAME(:,2,index_Max_1(index_Min_1(i))),'-o','LineWidth',2)
    string = append(['Combination ',num2str(index_Max_1(index_Min_1(i))),' Integral value: ',num2str(value_Max_1(index_Min_1(i))),'\n' ...
             'b[m]: ',num2str(filter_combination(index_Max_1(index_Min_1(i)),1)),' MAC[m]: ',num2str(filter_combination(index_Max_1(index_Min_1(i)),2)) ...
             ' h[m]: ',num2str(filter_combination(index_Max_1(index_Min_1(i)),3)),' s[m]: ',num2str(filter_combination(index_Max_1(index_Min_1(i)),4)),'\n' ...
             'C_{Lalpha}: ',num2str(C_L_alpha(index_Max_1(index_Min_1(i)))),'C_{L0}: ',num2str(C_L_0(index_Max_1(index_Min_1(i)))),'\n' ...
             'Reference surf [m^2]: ',num2str(surf_vec(index_Max_1(index_Min_1(i)))),' Normalized weight [m^2]: ',num2str(area_vec(index_Max_1(index_Min_1(i)))),'\n']);
    fprintf(string)
end
xlabel('\alpha','FontSize',14)
ylabel('C_L','FontSize',14)
xlim([-15 9.4])
set(gca,'FontSize',14)
grid on

figure(3)
plot(1./value_Max_1(1:count),area_vec(index_Max_1(1:count)),'o','MarkerFaceColor','#0072BD','LineWidth',2)
hold on
plot(1./value_Max_1(index_Min_1(1)),area_vec(index_Max_1(index_Min_1(1))),'o','MarkerFaceColor','#D95319','LineWidth',2)
plot(1./value_Max_1(index_Min_1(2)),area_vec(index_Max_1(index_Min_1(2))),'o','MarkerFaceColor','#EDB120','LineWidth',2)
plot(1./value_Max_1(index_Min_1(3)),area_vec(index_Max_1(index_Min_1(3))),'o','MarkerFaceColor','#7E2F8E','LineWidth',2)
grid on
xlabel('I^{-1}','FontSize',14)
ylabel('A [m^2]','FontSize',14)
% xlim([7.5 25])
% ylim([25 325])
set(gca,'FontSize',14)
legend('Combinations','Optimum 1','Optimum 2','Optimum 3','Location','Best')

function [C_L_alpha,C_L_0,valid] = liftCurve(data_APAME)
[a,b,c] = size(data_APAME);
for i = 1:c
    C_L_APAME = data_APAME(:,2,i);
    alpha_APAME = data_APAME(:,1,i);
    [~,R1]=fit(alpha_APAME,C_L_APAME,'poly1');
    [~,R2]=fit(alpha_APAME,C_L_APAME,'poly2');
    C_L_alpha(i) = (C_L_APAME(end)-C_L_APAME(1))/((alpha_APAME(end)-alpha_APAME(1))*pi/180);
    C_L_0(i) = interp1(alpha_APAME,C_L_APAME,0);
    if R1.rsquare > 0.9995 && C_L_alpha(i) < 15 && C_L_alpha(i) >= 0
        valid(i) = 1;
    else
        valid(i) = 0;
    end
end
end
function [value_Max,index_Max] = nMaxValues(variable,n)
    value_Max = zeros(n,1);
    index_Max = zeros(n,1);
    for j = 1:n
        [val,ind] = max(variable);
        variable(ind) = NaN;
        value_Max(j) = val;
        index_Max(j) = ind;
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
