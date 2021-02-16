%% set values for cantoni_LMI (Cantoni 2008)

nPool = 5;

% TODO: remove: data present in PoolModel (in secs and rad/sec!!!)
% % Pool model parameters
tau = [4, 2, 4, 4, 6]; % minutes
alpha = [6492, 2478, 6084, 5658, 7650]; % m^2
phi_wave = [0.48, 1.05, 0.48, 0.48, 0.42 ]; % rad/min (wave frequency)

% CVX suddenly decided to not be able to find its own function 'vec', so:
%addpath C:\Users\Jacob\Documents\MATLAB\CVX\cvx\functions\vec_
addpath ../functions_jacob/

% NOTE: The simulation was initially built for decentralized control using
% values from Cantoni 2007: Control of large scalle irrigation systems.
% That simulation still exists, but for distributed control I modified the
% simulation
% Values from Table 1 of the 2007 paper (for five pools):
% tau = [8,3,16,16,16]; % minutes; delay
% alpha = [22414, 11942, 43806, 43806, 43806]; % m^2
% phi_wave = [0.42, 0.74, 0.20, 0.20, 0.20]; % rad/min, dominant wave frequencies


% TODO: remove, zeta and omega_n also in PoolModel
% 3rd order model is used for simulation
% w_n = [0.65, 0.86, 0.448, 0.448, 0.448]; Wrong: computed wrongly
% We can just compute w_n now that we have a sound expression:

w_n = zeros(1,length(phi_wave));
for i=1:length(phi_wave)
    zeta(i) = 0.0151;  % no unit % From Gabriel's ini script
    w_n(i) = phi_wave(i)/sqrt(1-zeta(i)^2); % Literature survey sec. 3-3-1.
%     display(w_n(i));
end


% TODO: tune parameters / use code from simulation_lab
% Loop shaping weights parameters from [1]
kappa = [1.69, 6.47, 2.37, 2.21, 1.68];
phi = [113.64, 37.17, 86.96, 96.15, 113.64];
rho = [9.97, 3.26, 7.60, 8.47, 9.97];
% TODO: BB: eta seems unused, find out what it is for
eta = [130, 223, 183, 170, 153]; 

% %% Define the matrices