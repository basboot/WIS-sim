% resample for HIL

%% Parameters
clear all; 

%% load pool parameters for lab setup
lab_setup_values;

SPS = 8;



for i = 1:nPool

    % Now construct discrete-time versions of the dynamic control matrices
    h = 1/60 * 1/SPS; % Sampling time - STILL NEEDS TO BE PROPERLY SELECTED <------------
 
    ddelay(i) = round(tau(i)/h); % Discrete delay
    
    % first order
    P{i} = tf([1],[alpha(i) 0 ]);

    % Discretize the plant models
    Pd{i} = c2d(P{i}, h, 'zoh');

end

%% pool 0 %%
P0 = tf([1],[Wis.area0 0 ]);
Pd0 = c2d(P0, h, 'zoh');
