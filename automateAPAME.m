clear; close all; clc;
% APAME path
pathAPAME = '"F:\Universidad\Master_in_Aerospace_Engineering\TFM\Apame_v140915\ApameSolver\bin\apame_win64.exe"';
% Resolution study data
sweep_wing = [3 10 24 53 110];
sweep_winglet = [0 2 6 14 30];
sweep_airfoil = [29 29 29 29 29];
alpha_vec = [-2 4 12.3 14.5 16.5 18.5 21];
data_APAME = zeros(length(alpha_vec),4,length(sweep_airfoil));

for jindex = 1:length(sweep_airfoil)
factor = 0.0254;
h = 8*factor;                                    % [m] Gap between the wings
b = 40*factor;                                      % [m] Wingspan
R = 0.15*8*factor;                                    % [m] Cant radius
c_r1 = 8*factor;                                   % [m] Chord at the root of airfoil 1
c_r2 = 8*factor;                                   % [m] Chord at the root of airfoil 2
c_t1 = 8*factor;                                 % [m] Chord at the tip of airfoil 1
c_t2 = 8*factor;                                 % [m]Chord at the tip of airfoil 2
s = 0;                                      % [m] Stagger
N_bu = sweep_wing(jindex);                                   % Number of spanwise sections upper wing
N_bl = sweep_wing(jindex);                                   % Number of spanwise sections lower wing
N_bw = sweep_winglet(jindex);                                   % Number of spanwise sections winglet
Finite = 0;                                 % Finite TE
N_el = sweep_airfoil(jindex);
N = 5;
af_1r = '0012';
af_2r = '0012';
af_1t = '0012';
af_2t = '0012';
af_w = '0003';

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

aux_count = 0;
for indAoA = alpha_vec
    aux_count = aux_count +1;
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
    NoPan = BoxWingAPAME(Slices_Wing,h,name,airspeed,density,pressure,mach,cases,wingspan,MAC,surf,origin,method,err,colldepth,farfield,collcalc,velorder,results,requests);
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
    data_APAME(aux_count,1,jindex) = AoA; data_APAME(aux_count,2,jindex) = value_CL; data_APAME(aux_count,3,jindex) = value_CD; data_APAME(aux_count,4,jindex) = NoPan;
end
end
save('data_APAME_spanwise.mat','data_APAME')

%% Plots
clear; close all; clc;
load('data_APAME_spanwise.mat')
load('Validation_Data.mat')

figure()
errorbar(C_L_Exp(:,1),C_L_Exp(:,2),0.02*ones(length(C_L_Exp(:,1)),1),0.02*ones(length(C_L_Exp(:,1)),1),1*ones(length(C_L_Exp(:,1)),1),1*ones(length(C_L_Exp(:,1)),1),'.k','LineWidth',2,'MarkerSize',1);
hold on
plot(data_APAME(:,1,1),data_APAME(:,2,1),'o','color','#0072BD','LineWidth',2); grid on
plot(data_APAME(:,1,2),data_APAME(:,2,2),'o','color','#D95319','LineWidth',2); 
plot(data_APAME(:,1,3),data_APAME(:,2,3),'o','color','#EDB120','LineWidth',2);
plot(data_APAME(:,1,4),data_APAME(:,2,4),'o','color','#7E2F8E','LineWidth',2);
plot(data_APAME(:,1,5),data_APAME(:,2,5),'o','color','#77AC30','LineWidth',2);
xlabel('\alpha [\circ{}]','FontSize',14)
ylabel('$C_L$','FontSize',14,'Interpreter','latex')
set(gca,'FontSize',14)
legend('Experimental','1972 panels','3828 panels','7540 panels','15196 panels','30276 panels','Location','Best')
grid on

figure()
errorbar(C_L_Exp(1:9,2),C_D_Exp(:,2),0.0005*ones(length(C_D_Exp(:,1)),1),0.0005*ones(length(C_D_Exp(:,1)),1),0.02*ones(length(C_D_Exp(:,1)),1),0.02*ones(length(C_D_Exp(:,1)),1),'.k','LineWidth',2,'MarkerSize',1);
hold on
plot(data_APAME(:,2,1),data_APAME(:,3,1),'o','LineWidth',2);
plot(data_APAME(:,2,2),data_APAME(:,3,2),'o','LineWidth',2); 
plot(data_APAME(:,2,3),data_APAME(:,3,3),'o','LineWidth',2); 
plot(data_APAME(:,2,4),data_APAME(:,3,4),'o','LineWidth',2); 
plot(data_APAME(:,2,5),data_APAME(:,3,5),'o','LineWidth',2);
xlabel('$C_L$','FontSize',14,'Interpreter','latex')
ylabel('$C_D$','FontSize',14,'Interpreter','latex')
set(gca,'FontSize',14)
legend('Experimental','1972 panels','3828 panels','7540 panels','15196 panels','30276 panels','Location','Best')
grid on

err1 = abs((data_APAME(:,2,1)-C_L_Exp(1:8,2))./max(C_L_Exp(:,2))); av1 = mean(err1);
err2 = abs((data_APAME(:,2,2)-C_L_Exp(1:8,2))./max(C_L_Exp(:,2))); av2 = mean(err2);
err3 = abs((data_APAME(:,2,3)-C_L_Exp(1:8,2))./max(C_L_Exp(:,2))); av3 = mean(err3);
err4 = abs((data_APAME(:,2,4)-C_L_Exp(1:8,2))./max(C_L_Exp(:,2))); av4 = mean(err4);
err5 = abs((data_APAME(:,2,5)-C_L_Exp(1:8,2))./max(C_L_Exp(:,2))); av5 = mean(err5);

figure()
alpha = data_APAME(:,1,1);
plot(alpha,err1*100,'o','color','#0072BD','LineWidth',2)
hold on
plot(alpha,av1*100*ones(length(alpha),1),'--','color','#0072BD','LineWidth',2,'HandleVisibility','off')
plot(alpha,err2*100,'o','color','#D95319','LineWidth',2)
plot(alpha,av2*100*ones(length(alpha),1),'--','color','#D95319','LineWidth',2,'HandleVisibility','off')
plot(alpha,err3*100,'o','color','#EDB120','LineWidth',2)
plot(alpha,av3*100*ones(length(alpha),1),'--','color','#EDB120','LineWidth',2,'HandleVisibility','off')
plot(alpha,err4*100,'o','color','#7E2F8E','LineWidth',2)
plot(alpha,av4*100*ones(length(alpha),1),'--','color','#7E2F8E','LineWidth',2,'HandleVisibility','off')
plot(alpha,err5*100,'o','color','#77AC30','LineWidth',2)
plot(alpha,av5*100*ones(length(alpha),1),'--','color','#77AC30','LineWidth',2,'HandleVisibility','off')
xlabel('$\alpha [\circ]$','FontSize',14,'Interpreter','latex')
ylabel('Relative error [\%]','FontSize',14,'Interpreter','latex')
set(gca,'FontSize',14)
legend('1972 panels','3828 panels','7540 panels','15196 panels','30276 panels','Location','Best')
grid on