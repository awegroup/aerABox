% This function aims to read the data coming directly from OpenFOAM to
% check the convergence of the variables
%       Author  : Gabriel Buendia
%       Version : 1
% Inputs:
%       name -> String specifying the name (and path) of the file to be
%       read
%       numIt -> Number of iterations for the Cauchy residual
% Outputs:
%       forces_c_val -> Matrix containing the forces coeffients at
%                      each timestep (timestepx3)
%       res          -> Vector containing the Cauchy residuals of the
%                      aerodynamic forces (1x3)
function [forces_c_val,res] = ofScanLogSimpleFOAM(name,numIt)
fileID = fopen(name,'r');
string_full = fscanf(fileID,'%c');
fclose(fileID);
Ux_s = 'Ux, Initial residual = ';
Ux_index = strfind(string_full,Ux_s);
Uy_s = 'Uy, Initial residual = ';
Uy_index = strfind(string_full,Uy_s);
Uz_s = 'Uz, Initial residual = ';
Uz_index = strfind(string_full,Uz_s);
p_s = 'p, Initial residual = ';
p_index = strfind(string_full,p_s);
omega_s = 'omega, Initial residual = ';
omega_index = strfind(string_full,omega_s);
k_s = 'k, Initial residual =';
k_index = strfind(string_full,k_s);
forces_s = 'Sum of forces';
forces_index = strfind(string_full,forces_s)+21;
Ux_val = zeros(1,length(Ux_index));
Uy_val = zeros(1,length(Uy_index));
Uz_val = zeros(1,length(Uz_index));
p_val = zeros(2,length(p_index)/2);
omega_val = zeros(1,length(omega_index));
k_val = zeros(1,length(k_index));
forces_c_val = zeros(3,length(forces_index));
for i = 1:length(Ux_index)
    aux_ind = 1;
    while ~strcmp(string_full(Ux_index(i)+length(Ux_s)+aux_ind),char(44))
        aux_ind = aux_ind+1;
    end
    aux_ind = aux_ind-1;
    Ux_val(i) = str2num(string_full(Ux_index(i)+length(Ux_s):Ux_index(i)+length(Ux_s)+aux_ind));
end
for i = 1:length(Uy_index)
    aux_ind = 1;
    while ~strcmp(string_full(Uy_index(i)+length(Uy_s)+aux_ind),char(44))
        aux_ind = aux_ind+1;
    end
    aux_ind = aux_ind-1;
    Uy_val(i) = str2num(string_full(Uy_index(i)+length(Uy_s):Uy_index(i)+length(Uy_s)+aux_ind));
end
for i = 1:length(Uz_index)
    aux_ind = 1;
    while ~strcmp(string_full(Uz_index(i)+length(Uz_s)+aux_ind),char(44))
        aux_ind = aux_ind+1;
    end
    aux_ind = aux_ind-1;
    Uz_val(i) = str2num(string_full(Uz_index(i)+length(Uz_s):Uz_index(i)+length(Ux_s)+aux_ind));
end
index = 1;
for i = 1:length(p_index)
    aux_ind = 1;
    while ~strcmp(string_full(p_index(i)+length(p_s)+aux_ind),char(44))
        aux_ind = aux_ind+1;
    end
    aux_ind = aux_ind-1;
    if mod(i,2) == 1
        p_val(1,index) = str2num(string_full(p_index(i)+length(p_s):p_index(i)+length(Ux_s)+aux_ind));
    else
        p_val(2,index) = str2num(string_full(p_index(i)+length(p_s):p_index(i)+length(Ux_s)+aux_ind));
        index = index+1;
    end
end
for i = 1:length(omega_index)
    aux_ind = 1;
    while ~strcmp(string_full(omega_index(i)+length(omega_s)+aux_ind),char(44))
        aux_ind = aux_ind+1;
    end
    aux_ind = aux_ind-1;
    omega_val(i) = str2num(string_full(omega_index(i)+length(omega_s):omega_index(i)+length(omega_s)+aux_ind));
end
for i = 1:length(k_index)
    aux_ind = 1;
    while ~strcmp(string_full(k_index(i)+length(k_s)+aux_ind),char(44))
        aux_ind = aux_ind+1;
    end
    aux_ind = aux_ind-1;
    k_val(i) = str2num(string_full(k_index(i)+length(k_s):k_index(i)+length(k_s)+aux_ind));
end

for i = 1:length(forces_index)
    aux_ind = 1;
	for j = 1:3
		switch j 
			case 1
				while ~strcmp(string_full(forces_index(i)+length(forces_s)+aux_ind),char(32))
					aux_ind = aux_ind+1;
				end	
				aux_ind2 = aux_ind-1;
				aux_ind = aux_ind+1;
			case 2 
				while ~strcmp(string_full(forces_index(i)+length(forces_s)+aux_ind),char(32))
					aux_ind = aux_ind+1;
				end
				aux_ind3 = aux_ind-1;
				aux_ind = aux_ind + 1;
			otherwise
				while ~strcmp(string_full(forces_index(i)+length(forces_s)+aux_ind),char(41))
					aux_ind = aux_ind+1;
				end
				aux_ind4 = aux_ind-1;
        end
    end
    forces_c_val(1,i) = str2num(string_full(forces_index(i)+length(forces_s):forces_index(i)+length(forces_s)+aux_ind2));
	forces_c_val(2,i) = str2num(string_full(forces_index(i)+length(forces_s)+aux_ind2+1:forces_index(i)+length(forces_s)+aux_ind3));
	forces_c_val(3,i) = str2num(string_full(forces_index(i)+length(forces_s)+aux_ind3+1:forces_index(i)+length(forces_s)+aux_ind4));
end

figure()
semilogy(Ux_val,'Linewidth',2)
hold on
semilogy(Uy_val,'Linewidth',2)
semilogy(Uz_val,'Linewidth',2)
semilogy(p_val(2,:),'Linewidth',2)
semilogy(omega_val,'Linewidth',2)
semilogy(k_val,'Linewidth',2)
set(gca,'FontSize',14)
xlabel('Iterations','FontSize',14)
ylabel('Residuals','FontSize',14)
leg = legend('U_x','U_y','U_z','p','\omega','k','Location','Best');
leg.NumColumns = 2;
grid on

load("Flight_Conditions.mat")
figure()
semilogy(abs(forces_c_val(3,:))/Q,'Linewidth',2)
hold on
semilogy(abs(forces_c_val(2,:))/Q,'Linewidth',2)
semilogy(abs(forces_c_val(1,:))/Q,'Linewidth',2)
set(gca,'FontSize',14)
xlabel('Iterations','FontSize',14)
ylabel('|Aerodynamic coefficients|','FontSize',14)
legend('C_L','C_S','C_D','Location','Best')
grid on

forces_c_val = forces_c_val';
fprintf('Drag coefficient\n')
resDrag = cauchyResidual(forces_c_val(:,1)/Q,numIt);
fprintf('Sideforce coefficient\n')
resSide = cauchyResidual(forces_c_val(:,2)/Q,numIt);
fprintf('Lift coefficient\n')
resLift = cauchyResidual(forces_c_val(:,3)/Q,numIt);
res = [resDrag resSide resLift];
end

function res = cauchyResidual(forces,numIt)
counter = 0;
summatory = 0;
while counter <= numIt - 1
    diff = abs(forces(length(forces)-counter)-forces(length(forces)-counter-1));
    summatory = summatory + diff;
    counter = counter + 1;
end
res = summatory/(numIt-1);
fprintf('The Cauchy residual is %d\n',res)
end