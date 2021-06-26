%% set values for cantoni_LMI (lab setup)

%% Load identification results
try
    load('../identification/identification.mat');
catch
    assert(false, "File 'identification.mat' does not exist. Run identification first.");
end

% TODO: CHECK /seconds <=> /minutes

nPool = 3; % The number of pools. 

tuning_process = 1; % Set to 1 when in the process of tuning for bode plots
tuning_version = '1'; % used to generate different figure names for different setups

%% Loop shaping weights parameters from [1] % NEED CHANGES FOR THE LAB SETUP


% With correct values
tuning_version = 'valid';
kappa = [0.010, 0.03, 0.01]; % kappa_i is used to set the loop-gain bandwidth – this should also sit below (1/?i) rad/min, because of the delay which is not reflected in Li
phi = [64, 60, 35]; % phi_i is used to introduce phase lead in the cross-over region to reduce the roll-off rate for stability and robustness
rho = [1, 0.5, 0.01]; % rho_i is used to provide additional roll-off beyond the loop-gain bandwidth to ensure sufficiently low gain at the (un-modelled) dominant wave frequency

%2
% version used in HIL
tuning_version = 'exceeds_phi_wave';
kappa = [0.3, 0.5, 0.3]; % kappa_i is used to set the loop-gain bandwidth – this should also sit below (1/?i) rad/min, because of the delay which is not reflected in Li
phi = [10, 10, 10]; % phi_i is used to introduce phase lead in the cross-over region to reduce the roll-off rate for stability and robustness
rho = [0.1, 0.1, 0.1]; % rho_i is used to provide additional roll-off beyond the loop-gain bandwidth to ensure sufficiently low gain at the (un-modelled) dominant wave frequency


addpath ../functions_jacob/

for i = 1:nPool

    % Pool model parameters
    tau(i) = PoolModel(i).tau/60; % minutes 
    alpha(i) = PoolModel(i).alpha; % m^2 actual lab setup values
    phi_wave(i) = PoolModel(i).phi_wave*60; % rad/min (wave frequency) 
    zeta(i) = PoolModel(i).zeta;
    w_n(i) = PoolModel(i).omega_n;
    
%     fprintf("calculated using seconds: zeta %d omega_n %d\n", zeta(i), w_n(i));
%     zeta(i) = 0.0151;  % no unit % From Gabriel's ini script
%     w_n(i) = phi_wave(i)/sqrt(1-zeta(i)^2); % Literature survey sec. 3-3-1.
%     fprintf("calculated using minutes: zeta %d omega_n %d\n", zeta(i), w_n(i));

    W{i} = tf([kappa(i)*phi(i) kappa(i)], [rho(i) 1 0]); % shaping weight 

    % Shaping weight tuning based using the procedure described in [1]
    % tuning rules related to kappa(1)

    s = tf('s');
    L{i} = W{i}*(1/s*alpha(i)); % local loop-gain for W1
    freq(i) = (2*pi)/(tau(i)*60); % 1/tau(1) in rad/s


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
    
    bw_gain_phi(i) = evalfr(L{i}, phi_wave(i)/60); % magnitude/Gain, used to check if the bandwidth is < 1/tau(i)
    if bw_gain_phi(i) >= 0.7079 % 0.7079 corresponds to -3 dB in magnitude
        fprintf('The bandwidth of L%d is too large (bw_gain phi is: %d. \n Retune kappa before proceeding\n', i, bw_gain_phi(i));
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
        xline(phi_wave(i)/60, '--r');
        legend(sprintf('L%d', i),sprintf('1/s*alpha'), '1/\tau', '\phi_{wave}'); 
        grid on;
        saveFigureEps(sprintf('bode_weight_%d_version_%s',i,tuning_version));
        
        P{i} = 1/(alpha(i)*s);
        CLW{i} = feedback(P{i}*W{i},1);

        figure(); 
        step(CLW{i})
        title(sprintf('CLW %d', i)); 
        fprintf('Poles %d', i);
        pole(CLW{i})
        

    end % if in tuning process
end


%%%% CLEANUP

%kappa = [0.002, 0.003]; % kappa_i is used to set the loop-gain bandwidth – this should also sit below (1/?i) rad/min, because of the delay which is not reflected in Li
% BB: just a random changes to find values that work => TODO: look at this later!
% kappa = [0.01, 0.03, 0.01]; % kappa_i is used to set the loop-gain bandwidth – this should also sit below (1/?i) rad/min, because of the delay which is not reflected in Li
% phi = [45, 40, 35]; % phi_i is used to introduce phase lead in the cross-over region to reduce the roll-off rate for stability and robustness
% rho = [0.01, 0.01, 0.01]; % rho_i is used to provide additional roll-off beyond the loop-gain bandwidth to ensure sufficiently low gain at the (un-modelled) dominant wave frequency
% eta = [130, 223]; % Used for plot scaling in the same way as done in [1]


% % retune
% kappa = [0.00018292, 0.00069301, 4.4902e-05]; % kappa_i is used to set the loop-gain bandwidth – this should also sit below (1/?i) rad/min, because of the delay which is not reflected in Li
% phi = [1.6e+02, 53, 4.4e+02]; % phi_i is used to introduce phase lead in the cross-over region to reduce the roll-off rate for stability and robustness
% rho = [0.32, 0.16, 0.42]; % rho_i is used to provide additional roll-off beyond the loop-gain bandwidth to ensure sufficiently low gain at the (un-modelled) dominant wave frequency


% % tune slower
% kappa = [0.00001, 0.00001, 0.00001]; % kappa_i is used to set the loop-gain bandwidth – this should also sit below (1/?i) rad/min, because of the delay which is not reflected in Li
% phi = [5000, 5000, 5000]; % phi_i is used to introduce phase lead in the cross-over region to reduce the roll-off rate for stability and robustness
% rho = [50, 50, 50]; % rho_i is used to provide additional roll-off beyond the loop-gain bandwidth to ensure sufficiently low gain at the (un-modelled) dominant wave frequency

% BB: retune to make it slower
% kappa = [0.01, 0.03, 0.01]; % kappa_i is used to set the loop-gain bandwidth – this should also sit below (1/?i) rad/min, because of the delay which is not reflected in Li
% phi = [45, 40, 35]; % phi_i is used to introduce phase lead in the cross-over region to reduce the roll-off rate for stability and robustness
% rho = [0.01, 0.01, 0.01]; % rho_i is used to provide additional roll-off beyond the loop-gain bandwidth to ensure sufficiently low gain at the (un-modelled) dominant wave frequency

% % BB: retune more simular to cantoni, ignoring phi_wave
% version sent to Manuel and Gabriel
% kappa = [0.2, 0.3, 0.2]; % kappa_i is used to set the loop-gain bandwidth – this should also sit below (1/?i) rad/min, because of the delay which is not reflected in Li
% phi = [10, 10, 10]; % phi_i is used to introduce phase lead in the cross-over region to reduce the roll-off rate for stability and robustness
% rho = [0.1, 0.1, 0.1]; % rho_i is used to provide additional roll-off beyond the loop-gain bandwidth to ensure sufficiently low gain at the (un-modelled) dominant wave frequency


% BB: retune using auti PI (rho added very small)
% tuning_version = '3';
% kappa = [0.65, 0.45, 0.85]; % kappa_i is used to set the loop-gain bandwidth – this should also sit below (1/?i) rad/min, because of the delay which is not reflected in Li
% phi = [0.75, 0.65, 0.65]; % phi_i is used to introduce phase lead in the cross-over region to reduce the roll-off rate for stability and robustness
% rho = [0.0001, 0.0001, 0.0001]; % rho_i is used to provide additional roll-off beyond the loop-gain bandwidth to ensure sufficiently low gain at the (un-modelled) dominant wave frequency
% 