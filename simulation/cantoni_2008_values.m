%% set values for cantoni_LMI (Cantoni 2008)

nPool = 5; % added 6th pool to test nPool script

% TODO: remove: data present in PoolModel (in secs and rad/sec!!!)
% % Pool model parameters
tau = [4, 2, 4, 4, 6, 4]; % minutes
alpha = [6492, 2478, 6084, 5658, 7650, 6492]; % m^2
phi_wave = [0.48, 1.05, 0.48, 0.48, 0.42 , 0.48]; % rad/min (wave frequency)

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
kappa = [1.69, 6.47, 2.37, 2.21, 1.68, 1.69];
phi = [113.64, 37.17, 86.96, 96.15, 113.64, 113.64];
rho = [9.97, 3.26, 7.60, 8.47, 9.97, 9.97];
% BB: eta is only used for scaling of the plots
eta = [130, 223, 183, 170, 153, 130]; 

% %% Define the matrices

%% check Cantoni bode plots
tuning_process = 1;

for i = 1:1

    W{i} = tf([kappa(i)*phi(i) kappa(i)], [rho(i) 1 0]); % shaping weight 

    % Shaping weight tuning based using the procedure described in [1]
    % tuning rules related to kappa(1)

    s = tf('s');
    L{i} = W{i}*(1/s*alpha(i)); % local loop-gain for W1
    freq(i) = (2*pi)/(tau(i)); % 1/tau(1) in rad/min


    % Tuning phi and rho:
    disp('phi is tuned to introduce phase lead in the cross-over region to reduce the roll-off rate for stability and robustness');
    disp('rho is tuned to provide additional roll-off beyond the loop-gain bandwidth to ensure sufficiently low gain at the (un-modelled) dominant wave frequency');

    % Tuning kappa
    % bw_gain should be less than 0.7079 at the 'freq', then the bandwidth is < 1/tau(i)
    bw_gain_tau(i) = evalfr(L{i}, freq(i)); % magnitude/Gain, used to check if the bandwidth is < 1/tau(i)
    if bw_gain_tau(i) >= 0.7079 % 0.7079 corresponds to -3 dB in magnitude
        fprintf('The bandwidth of L%d is too large (bw_gain tau is: %d. \n Retune kappa before proceeding\n', i, bw_gain_tau(i));
        %return;
    else
        fprintf('The bandwidth of L%d is good (bw_gain tau is %d \n', i, bw_gain_tau(i));
    end
    
    bw_gain_phi(i) = evalfr(L{i}, phi_wave(i)); % magnitude/Gain, used to check if the bandwidth is < 1/tau(i)
    if bw_gain_phi(i) >= 0.7079 % 0.7079 corresponds to -3 dB in magnitude
        fprintf('The bandwidth of L%d is too large (bw_gain phi is: %d. \n Retune kappa(1) before proceeding\n', i, bw_gain_phi(i));
        %return;
    else
        fprintf('The bandwidth of L%d is good (bw_gain phi is %d \n', i, bw_gain_phi(i));
    end

    if tuning_process  == 1
        % Make bode plots for L1 and L2
        figure(); 
        bode(L{i}); % bodeJL(W1,'Plantname');
        hold on; 
        bode(1/s*alpha(i)); 
        xline(freq(i),'--b');
        xline(phi_wave(i), '--r');
        legend(sprintf('L%d', i),sprintf('1/s*alpha'), '1/\tau', '\phi_{wave}'); 
        grid on;
        
        P{i} = 1/(alpha(i)*s);
        CLW{i} = feedback(P{i}*W{i},1);

        figure(); 
        step(CLW{i})
        title(sprintf('CLW %d', i)); 
        fprintf('Poles %d', i);
        pole(CLW{i})

    end % if in tuning process
end