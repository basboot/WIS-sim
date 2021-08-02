% Author: Jacob Lont
% Date: 10-2019
% Initialization script for the simulink file cantonisim_distributed.slx
% Update: Not used anymore at this point. Just run cantoni_LMI.m once
% before using the simulation.

%% DISTRIBUTED CONTORL

% M-file used in combination with my simulink file in order to recreate a
% water irrigation system simulation based on Cantoni et al. (2007):
% Control of Large-Scale Irrigation Networks.

clc; clear all;

% Values from Table 1 of the paper (for five pools):
tau = [8,3,16,16,16]; % minutes; delay
alpha = [22414, 11942, 43806, 43806, 43806]; % m^2
phi_wave = [0.42, 0.74, 0.20, 0.20, 0.20]; % rad/min, dominant wave frequencies


%% Constructing the controller
% Values NOT given in the paper (only for identical pools):
kappa = [7.72, 0, 0, 0, 0];
phi = [128 0 0 0 0];
rho = [15.2 0 0 0 0];

% 1st order model is used for controller design
% s = tf('s');
% C1 = kappa(1)*(1+ s*phi(1))/ (s*(1+s*rho(1)));
% bode(C1)

%% 3rd order model is used for simulation
% Try to compute the values (2 unknowns)
w_n = [0.65, 0.86, 0.448, 0.448, 0.448];
zeta = 0.0151;  % no unit % From Gabriel's ini script

%for i=1:5
%phi_wave_computed(i) = [w_n(i)^2 * (1-zeta^2)];
%end
% display(phi_wave_computed);

% Construct the models:


%% Create generalized plants
s = tf('s');

for i=1:5

    G11 = 0;
    G12 = [0 0 1];
    G13 = 1;

    G21 = [1/(s*alpha(i)); 0];
    G22 = [1 1/(s*alpha(i)) -exp(-s*tau(i))/(s*alpha(i));
           0 0 0];
	G23 = [-exp(-s*tau(i))/(s*alpha(i)); 1];
    
    G31 = 1/(s*alpha(i));
    G32 = [1 1/(s*alpha(i)) -exp(-s*tau(i))/(s*alpha(i))];
    G33 = [-exp(-s*tau(i))/(s*alpha(i))];
    
    Gtemp = [G11 G12 G13;
             G21 G22 G23;
             G31 G32 G33];

    switch i
        case 1
            G1 = Gtemp;
        case 2
            G2 = Gtemp;
        case 3
            G3 = Gtemp;
        case 4
            G4 = Gtemp;
        case 5
            G5 = Gtemp;
    end
end

% Clean up workspace
clear G11 G12 G13 G21 G22 G23 G31 G32 G33 Gtemp;


%% Compute a stabilizing H?-optimal controller K for the plant P.

G = [G1 G2 G3 G4 G5]; % Cell array of the models Gi(s)

% Create a matrix to store the decentralized controllers
C(1) = kappa(1)*(1+ s*phi(1))/ (s*(1+s*rho(1))); % USe the same controller for all pools

% Create a combined weight matrix
Wm_mini = [eye(3) zeros(3,2); 
      zeros(2,3) C(1)*eye(2)];% Weight matrix for one plant
Wm = blkdiag(Wm_mini, Wm_mini, Wm_mini, Wm_mini, Wm_mini); % 25x25 

% Weighted generalized plant
Gw = G*Wm;

% padeGw = pade(Gw,3);
padeG1 = ss(pade(G1,1)*Wm_mini); % 1st order pade approximation
size(padeG1)

% H_inf synthesis
nmeas = 1;
ncont = 2;
% [K,CL,gamma] = hinfsyn(padeGw,nmeas,ncont);
[K,CL,gamma] = hinfsyn(padeG1,nmeas,ncont)

% https://nl.mathworks.com/help/robust/ref/hinfsyn.html
% nmeas and ncont are the number of signals in y and u, respectively. 
% y and u are the last outputs and inputs of P, respectively. 
% hinfsyn returns a controller K that stabilizes P and has the same number of states. 
% The closed-loop system CL = lft(P,K) achieves the performance level gamma, which is the H? norm of CL

























