% This function computes the number of elements given an element size, some
% geometrical parameters and the aspect ratio of the desired cells
%       Author  : Gabriel Buendia
%       Version : 1
% Inputs:
%       el_size  -> Size of the element [m]
%       b        -> Wingspan [m]
%       R        -> Cant radius [m]
%       h        -> Height between wings [m]
%       s        -> Wing stagger [m]
%       c_r1     -> Chord at the root of the lower wing [m]
%       c_t1     -> Chord at the tip of the lower wing [m]
%       c_r2     -> Chord at the root of the upper wing [m]
%       c_t2     -> Chord at the tip of the upper wing [m]
%       N        -> Number of sections to discretize the cant region
%       N_bu     -> Number of spanwise sections upper wing
%       N_bl     -> Number of spanwise sections lower wing
%       N_bw     -> Number of spanwise sections winglet
%       AR       -> Aspect ratio of the desired cells
% Outputs:
%       N_LEls   -> Number of elements in the lower wing leading edge
%       N_LEus   -> Number of elements in the upper wing leading edge
%       N_AF     -> Number of elements present in the airfoil
%       N_R      -> Number of elements in the leading edge of the cant region
%       N_WL     -> Number of elements in the winglet
function [N_LEls, N_LEus, N_AF, N_R, N_WL] = computeNumberOfElements(el_size,b,R,h,s,c_r1,c_t1,c_r2,c_t2,N,N_bl,N_bu,N_bw,AR)
c_min = min([c_r1,c_r2,c_t1,c_t2]);
N_LEls = 1+ceil((b/2-R)/((N_bl+1)*el_size)/AR);
N_LEus = 1+ceil((b/2-R)/((N_bu+1)*el_size)/AR);
N_AF = 1+ceil(c_min/el_size);
N_R = 1+ceil(pi*R/(2*(N-1))/el_size/AR);
N_WL = 1+ceil((sqrt(h^2+s^2)-2*R)/((N_bw+1)*el_size)/AR);
end

