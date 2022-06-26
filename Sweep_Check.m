clear; close all; clc;
span_initial = 42.4694;
chord_initial = 3.5425;
surf_lim = 150.4478;
area_lim = 323.3616;
taper = 1;
fore_area = 0.5;
% sweep_b = linspace(span_initial/4,span_initial,10);
% sweep_c = linspace(chord_initial/4,chord_initial,10);
sweep_b = linspace(span_initial*4/10,span_initial*6/10,11);
sweep_c = linspace(chord_initial*2/10,chord_initial*3/10,11);
% sweep_h = linspace(0.2*span_initial/4,span_initial,10);
% sweep_s = linspace(0,chord_initial,10);
sweep_h = linspace(0.05,0.15,11);
sweep_s = linspace(-1,1,11);
valid = zeros(length(sweep_b),length(sweep_c),length(sweep_h),length(sweep_s));
combination = zeros(5000,4);
counter = 0;
for i = 1:length(sweep_b)
    for j = 1:length(sweep_c)
        for k = 1:length(sweep_h)
            for l = 1:length(sweep_s)
                b = sweep_b(i); wingspan = b;
                MAC = sweep_c(j); c_r1 = MAC; c_r2 = MAC; c_t1 = MAC; c_t2 = MAC;
                surf = 2*MAC*b;
%                 h = sweep_h(k);
%                 s = sweep_s(l);
                h = sweep_h(k)*b;
                s = sweep_s(l)*MAC;
                R = 0.15*sqrt(h^2+s^2);  
                N_bu = 3; N_bl = 3; N_bw = 0; N_el = 30; N = 5;                        
                Finite = 0;                               
                af_1r = '23012'; af_2r = '23012'; af_1t = '23012'; af_2t = '23012'; af_w = '0030';
                [Slices_Wing, Points_Wing] = ParametrizationBWFiniteTE(h,b,R,c_r1,c_r2,c_t1,c_t2,s,N_bu,N_bl,N_bw,N_el,N,Finite,af_1r,af_2r,af_1t,af_2t,af_w);
                name = 'A';
                AoA = 0; AoS = 0; cases = [AoA; AoS];
                height = 0;
                [T,~,p,~] = atmosisa(height);
                [rho,mu] = AirProperties(T-273.85, p/100, [],'rho', 'mu');
                gamma = 1.4; R_c = 287.04; density = rho; pressure = p; 
                airspeed = 60; Re = rho*MAC*airspeed/mu; mach = airspeed/sqrt(gamma*R_c*T);
                origin = [0 0 0]; method = 0; err = 1e-7; colldepth = 1e-7; farfield = 5; collcalc = 1; velorder = 2; results = 1; requests = zeros(1,13); requests(1) = 1;
                [NoPan,area] = BoxWingAPAME(Slices_Wing,h,name,airspeed,density,pressure,mach,cases,wingspan,MAC,surf,origin,method,err,colldepth,farfield,collcalc,velorder,results,requests);
                if (surf <= surf_lim) && (area <= area_lim) 
                    counter = counter+1;
                    valid(i,j,k,l) = 1;
                    combination(counter,1) = b;
                    combination(counter,2) = MAC;
                    combination(counter,3) = h;
                    combination(counter,4) = s;
                end
            end
        end
    end
end
save('SweepBounded2.mat') % 1 Is the initial sweep (commented)

%% Plots
clear; close all; clc;
load('SweepBounded2.mat') % 1 Is the initial sweep (commented)
figure()
histogram(combination(:,3)./combination(:,1))
xlabel('h/b','FontSize',14)
set(gca,'FontSize',14)
figure()
histogram(combination(:,4)./combination(:,2))
xlabel('s/c','FontSize',14)
set(gca,'FontSize',14)