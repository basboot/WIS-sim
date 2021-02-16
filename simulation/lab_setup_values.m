%% set values for cantoni_LMI (lab setup)

%% Load identification results
try
    load('../identification/identification.mat');
catch
    assert(false, "File 'identification.mat' does not exist. Run identification first.");
end

% TODO: CHECK /seconds <=> /minutes

nPool = 3; % The number of pools. 

%% Loop shaping weights parameters from [1] % NEED CHANGES FOR THE LAB SETUP
%kappa = [0.002, 0.003]; % kappa_i is used to set the loop-gain bandwidth – this should also sit below (1/?i) rad/min, because of the delay which is not reflected in Li
% BB: just a random changes to find values that work => TODO: look at this later!
kappa = [0.02, 0.03, 0.04]; % kappa_i is used to set the loop-gain bandwidth – this should also sit below (1/?i) rad/min, because of the delay which is not reflected in Li
phi = [45, 40, 35]; % phi_i is used to introduce phase lead in the cross-over region to reduce the roll-off rate for stability and robustness
rho = [0.01, 0.01, 0.01]; % rho_i is used to provide additional roll-off beyond the loop-gain bandwidth to ensure sufficiently low gain at the (un-modelled) dominant wave frequency
% eta = [130, 223]; % Used for plot scaling in the same way as done in [1]

tuning_process = 1; % Set to 1 when in the process of tuning for bode plots

addpath ../functions_jacob/

for i = 1:nPool

    % Pool model parameters
    tau(i) = PoolModel(i).tau/60; % minutes 
    alpha(i) = PoolModel(i).alpha; % m^2 actual lab setup values
    phi_wave(i) = PoolModel(i).phi_wave*60; % rad/min (wave frequency) 
    zeta(i) = PoolModel(i).zeta;
    w_n(i) = PoolModel(i).omega_n;
    

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
    bw_gain(i) = evalfr(L{i}, freq(i)); % magnitude/Gain, used to check if the bandwidth is < 1/tau(i)
    if bw_gain(i) >= 0.7079 % 0.7079 corresponds to -3 dB in magnitude
        fprintf('The bandwidth of L%d is too large (bw_gain(1) is: %d. \n Retune kappa(1) before proceeding\n', i, bw_gain(i));
        return;
    else
        fprintf('The bandwidth of L%d is good (bw_gain(1) is %d \n', i, bw_gain(i));
    end

    if tuning_process  == 1
        % Make bode plots for L1 and L2
        figure(); 
        bode(L{i}); % bodeJL(W1,'Plantname');
        hold on; 
        bode(1/s*alpha(i)); 
        xline(freq(i),'--r')
        legend(sprintf('L%d', i),sprintf('1/s*alpha(%d)', i), 'max bw'); 
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