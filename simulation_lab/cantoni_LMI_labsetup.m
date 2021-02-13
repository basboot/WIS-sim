%% Construction of a set of distributed controllers for a string of pools
% Source: Distributed controller design for open water channels (2008) [1]
% Yuping Li & Michael Cantoni

%% This file is for control design for the LAB SETUP
% Author: Jacob Lont with help from Gabriel A. Gleizer
% Date: 12-11-2019
% Last modified: 09-07-2020
% Status: The file is modified for the lab setup case, although parameters
% still need to be updated and next simulations tests have to be done.

% Todo: 
% - Measure the time constants of time plant using a step input
% - Check the sensitivity and complementary sensitivity functions for
% shaping weight tuning. CVX cannot solve the problem at this point.


%% Parameters
clear all; % for debug purposes, can be removed later
% Pool model parameters
tau = [0.16, 0.1]; % minutes % NEEDS TO BE UPDATED, VALUES ARE GUESSED ON INSTINCT
alpha = [0.1853, 0.1187]; % m^2 actual lab setup values
phi_wave = [3, 3]; % rad/min (wave frequency) % NEEDS TO BE UPDATED; NOT THE ACTUAL VALUES

% CVX suddenly decided to not be able to find its own function 'vec', so:
addpath C:\Users\Jacob\Documents\MATLAB\CVX\cvx\functions\vec_

% 3rd order model is used for simulation, which needs a frequency:
% We can just compute w_n now that we have a sound expression:
zeta = 0.0151;  % no unit % From Gabriel's ini script
w_n = zeros(1,length(phi_wave));
for i=1:length(phi_wave)
    w_n(i) = phi_wave(i)/sqrt(1-zeta^2); % Literature survey sec. 3-3-1.
%     display(w_n(i));
end

%% Loop shaping weights parameters from [1] % NEED CHANGES FOR THE LAB SETUP
kappa = [0.002, 0.003]; % kappa_i is used to set the loop-gain bandwidth – this should also sit below (1/?i) rad/min, because of the delay which is not reflected in Li
phi = [45, 40]; % phi_i is used to introduce phase lead in the cross-over region to reduce the roll-off rate for stability and robustness
rho = [1, 1]; % rho_i is used to provide additional roll-off beyond the loop-gain bandwidth to ensure sufficiently low gain at the (un-modelled) dominant wave frequency
% eta = [130, 223]; % Used for plot scaling in the same way as done in [1]


% Continuous-time shaping weight construction and tuning
W1 = tf([kappa(1)*phi(1) kappa(1)], [rho(1) 1 0]); % shaping weight 1
W2 = tf([kappa(2)*phi(2) kappa(2)], [rho(2) 1 0]); % shaping weight 2


% Shaping weight tuning based using the procedure described in [1]
% tuning rules related to kappa(1)


tuning_process = 1; % Set to 1 when in the process of tuning for bode plots

clc;
s = tf('s');
L1 = W1*(1/s*alpha(1)); % local loop-gain for W1
L2 = W2*(1/s*alpha(2)); % local loop-gain for W2
freq(1) = (2*pi)/(tau(1)*60); % 1/tau(1) in rad/s
freq(2) = (2*pi)/(tau(2)*60); % 1/tau(2) in rad/s

% Tuning phi and rho:
disp('phi is tuned to introduce phase lead in the cross-over region to reduce the roll-off rate for stability and robustness');
disp('rho is tuned to provide additional roll-off beyond the loop-gain bandwidth to ensure sufficiently low gain at the (un-modelled) dominant wave frequency');

% Tuning kappa(1)
% bw_gain should be less than 0.7079 at the 'freq', then the bandwidth is < 1/tau(i)
bw_gain(1) = evalfr(L1, freq(1)); % magnitude/Gain, used to check if the bandwidth is < 1/tau(i)
if bw_gain(1) >= 0.7079 % 0.7079 corresponds to -3 dB in magnitude
    fprintf('The bandwidth of L1 is too large (bw_gain(1) is: %d. \n Retune kappa(1) before proceeding\n', bw_gain(1));
    return;
else
    fprintf('The bandwidth of L1 is good (bw_gain(1) is %d \n', bw_gain(1));
end

% Tuning kappa(2)
% bw_gain should be less than 0.7079 at the 'freq', then the bandwidth is < 1/tau(i)
bw_gain(2) = evalfr(L2, freq(2)); % magnitude/Gain, used to check if the bandwidth is < 1/tau(i)
if bw_gain(2) >= 0.7079
    fprintf('The bandwidth of L2 is too large (bw_gain(2) is: %d. \n Retune kappa(2) before proceeding\n', bw_gain(2));
    return;
else
    fprintf('The bandwidth of L2 is good (bw_gain(2) is %d \n', bw_gain(2));
end
if tuning_process  == 1
    % Make bode plots for L1 and L2
    figure(1); clf; bode(L1); % bodeJL(W1,'Plantname');
    hold on; bode(1/s*alpha(1)); legend('L1','1/s*alpha(1)'); grid on;

    figure(2); clf; bode(L2);
    hold on; bode(1/s*alpha(2)); legend('L2','1/s*alpha(2)'); grid on;

    % Reminder
    fprintf('Try to compare the bode plots with the ones for the Cantoni values for comparison of loop shapes.\n');




    % Check some stability properties
    % The shaping weight can be used a decentralized compensator, shown in [1]
    % So, one could state that the feedback loop of Wi, Pi should be stable
    P1 = 1/(alpha(1)*s);
    P2 = 1/(alpha(2)*s);
    CLW1 = feedback(P1*W1,1);
    CLW2 = feedback(P2*W2,1);

    % figure(3); clf; bode(CLWi)
    figure(4); clf; step(CLW1)
    pole(CLW1)

    figure(5); clf; step(CLW2)
    pole(CLW2)

end % if in tuning process
%% Define the matrices

N = 2; % The number of pools. 
for i=1:N
    % Create the pool specific matrices for pool i
    % Corresponding to eq. (4) of the paper.
    Atti = [0, 1/alpha(i), -1/alpha(i), 0;
            0, -2/tau(i), 4/tau(i), 0;
            0, 0, 0, 1;
            0, 0, 0, -1/rho(i) ];
    Atsi = [-1/alpha(i); 0; 0; 0];
    Btni = [0, -1/alpha(i), 0;
            0, 0, 0;
            0, 0, kappa(i)*phi(i)/rho(i);
            0, 0, kappa(i)*(rho(i)-phi(i))/rho(i)^2];
    Btui = [0; 0; kappa(i)*phi(i)/rho(i); kappa(i)*(rho(i)-phi(i))/rho(i)^2];

    Asti = [0, 0, 1, 0];
    Assi = [0];
    Bsni = [0, 0, 0];
    Bsui = [0];
    Ctzi = [-1 , 0, 0, 0;
             0, 0, 0, 0];
    Cszi = [0; 0];
    Ctyi = [-1, 0, 0, 0];
    Csyi = [0];

    Dzni = [1, 0, 0;
            0, 0, 0];
    Dzui = [0; 1];
    Dyni = [1, 0, 0];
    Dyui = [0];

    % Construct NiX and NiY
    CCD = [Ctzi, Cszi, Dzni];
    BBD = [Btui', Bsui', Dzui'];

    NiX = null(CCD);
    NiY = null(BBD);

    % Construct \Pi_i^X
    PiXi = [eye(4),      zeros(4,1),   zeros(4,3);
            Atti,        Atsi,         Btni;
            Asti,        Assi,         Bsni;
            zeros(1,4),  1,            zeros(1,3);
            Ctzi,        Cszi,         Dzni;
            zeros(3,4),  zeros(3,1),   eye(3)]       * NiX;

    % Construct \Pi_i^Y
    PiYi = [Atti',       Asti',        Ctzi';
            -eye(4),     zeros(4,1),   zeros(4,2);
            zeros(1,4),  -eye(1),      zeros(1,2);       % Comment this line only for pool 1
            Atsi',       Assi',        Cszi';
            zeros(2,4),  zeros(2,1),   -eye(2);
            Btni',       Bsni',        Dzni']       * NiY;

    switch i
        case 1          % Pool 1 specific matrices
            N1X = NiX;
            % Modify N1Y to comply with the empty (5th) column of PiY1
            N1Y = [NiY(1:4,:);  NiY(6:7,:)]; % Only for pool 1

            Att1 = Atti; Ats1 = Atsi;
            Btn1 = Btni; Btu1 = Btui; Bsu1 = Bsui;

            Ast1 = zeros(0,4); % Only for pool 1
            Ass1 = zeros(0,1); % Only for pool 1
            Bsn1 = zeros(0,3); % Only for pool 1
%             Bsu1 = zeros(0,1); % Only for pool 1  % Initial mistake?

            Ctz1 = Ctzi; Csz1 = Cszi; Cty1 = Ctyi; Csy1 = Csyi;
            Dzn1 = Dzni; Dzu1 = Dzui; Dyn1 = Dyni; Dyu1 = Dyui;

            % PiX1, PiY1 computation differs from rest of the pools
            PiX1 = [eye(4),     zeros(4,1),  zeros(4,3);
                    Att1,       Ats1,        Btn1;
                    Ast1,       Ass1,        Bsn1;
                    zeros(1,4), 1,           zeros(1,3);
                    Ctz1,       Csz1,        Dzn1;
                    zeros(3,4), zeros(3,1),  eye(3)]      * N1X;
            
            PiY1 = [Att1',      Ast1',        Ctz1';
                    -eye(4),    zeros(4,0),   zeros(4,2);
                    %0,  -zeros(1,0),  0;       % Only comment for pool 1
                    Ats1',      Ass1',        Csz1';
                    zeros(2,4), zeros(2,0),   -eye(2);
                    Btn1',      Bsn1',        Dzn1']      * N1Y;

%   The code used for these cases is still included in the other directory
%         case 2          % Pool 2 specific matrices
%         case 3          % Pool 3 specific matrices
%         case 4          % Pool 4 specific matrices
%         case 5          % Pool 5 (end-pool) specific matrices

        case N          % Pool 2 (end-pool) specific matrices
%             N5X = NiX;
            N2Y = NiY;
            % Modify N5X to comply with the empty (5th) column of PiX5
            % Remove row 5 of N5X
            N2X = [NiX(1:4,:);  NiX(6:8,:)]; % Only for pool 5

            Att2 = Atti; Ast2 = Asti; %Ass5 = Assi;
            Btn2 = Btni; Btu2 = Btui; Bsn2 = Bsni; Bsu2 = Bsui;

            Ats2 = zeros(4,0); % Only for pool N
            Ass2 = zeros(1,0); % Only for pool N
%             Ass5 = []; % This is what we use later on to Create S5, so
%             maybe we should also implement this here?
            Csz2 = zeros(2,0); % Only for pool N
            
            Ctz2 = Ctzi; Cty2 = Ctyi; Csy2 = Csyi;
            Dzn2 = Dzni; Dzu2 = Dzui; Dyn2 = Dyni; Dyu2 = Dyui;

            % PiXN, PiYN computation differs from rest of the pools
            PiX2 = [eye(4),     zeros(4,0),  zeros(4,3);
                    Att2,       Ats2,        Btn2;
                    Ast2,       Ass2,        Bsn2;
                    %zeros(1,4), 1,           zeros(1,3); % Only pool 5
                    Ctz2,       Csz2,        Dzn2;
                    zeros(3,4), zeros(3,0),  eye(3)]      * N2X;
            
            PiY2 = [Att2',      Ast2',        Ctz2';
                    -eye(4),    zeros(4,1),   zeros(4,2);
                    zeros(1,4)  -eye(1),      zeros(1,2);
                    Ats2',      Ass2',        Csz2';
                    zeros(2,4), zeros(2,1),   -eye(2);
                    Btn2',      Bsn2',        Dzn2']      * N2Y;

    end % switch case
end %for


%% Define the optimization problem for the string of pools

% gamma_sqr = 8.2944; % <= gamma = 2.88 is optimum from the paper
% Smallest that can be solved: 
% gamma_sqr = 16; % -> gamma = 4
% gamma_sqr = 15.4; % -> gamma = 3.9243 % Min value for the Cantoni values
%gamma_sqr = 8; % -> gamma = 2.8284 % 2 pools with cantoni parammeters
gamma_sqr = 200; % -> gamma = ..
disp(strcat('solving for \gamma = ', num2str(sqrt(gamma_sqr))));

% Use CVX to perform the optimization
cvx_begin sdp  % semi-definite programming
    variable Xtt1(4,4) semidefinite; % implies pos def
    variable Xtt2(4,4) semidefinite;
    variable Ytt1(4,4) semidefinite;
    variable Ytt2(4,4) semidefinite;
    variable X21(1,1) symmetric; % symmetric, but not pos def
    variable Y21(1,1) symmetric;
%     variable gamma_sqr nonnegative


% The problem is said to be convex in the paper, but it is not in practice.
% Therefor we minimize for 0, instead of gamma_sqr
    minimize 0  
    
%     minimize gamma_sqr
    subject to
        % Did not include: Xitt > 0, Yitt > 0, quite sure this is already
        % implicitly declared.
        Y21 <= 0; %Y32 <= 0; Y43 <= 0; Y54 <= 0; % < 0
        X21 <= 0; %X32 <= 0; X43 <= 0; X54 <= 0; % < 0
        % 
        [Xtt1, eye(4); eye(4), Ytt1] >= 0;
        [Xtt2, eye(4); eye(4), Ytt2] >= 0;
        % 
        [X21, -1; -1, Y21] <= 0;
        % 
        PiX1' * ...
            [zeros(4)    Xtt1        zeros(4,1)   zeros(4,2)  zeros(4,3);
             Xtt1        zeros(4)    zeros(4,1)   zeros(4,2)  zeros(4,3);
             zeros(1,4)  zeros(1,4)  X21          zeros(1,2)  zeros(1,3);
             zeros(2,4)  zeros(2,4)  zeros(2,1)   eye(2)      zeros(2,3);
             zeros(3,4)  zeros(3,4)  zeros(3,1)   zeros(3,2)  -gamma_sqr*eye(3)] * ...
             PiX1 <= 0;
        PiY1' * ...
            [zeros(4)    Ytt1        zeros(4,1)   zeros(4,2)  zeros(4,3);
             Ytt1        zeros(4)    zeros(4,1)   zeros(4,2)  zeros(4,3);
             zeros(1,4)  zeros(1,4)  Y21          zeros(1,2)  zeros(1,3);
             zeros(2,4)  zeros(2,4)  zeros(2,1)   eye(2)      zeros(2,3);
             zeros(3,4)  zeros(3,4)  zeros(3,1)   zeros(3,2) -(1/gamma_sqr)*eye(3)] * ...
             PiY1 >= 0;
         % 2 (different, Nth pool)
         PiX2' * ...
            [zeros(4)    Xtt2        zeros(4,1)    zeros(4,2)  zeros(4,3);
             Xtt2        zeros(4)    zeros(4,1)    zeros(4,2)  zeros(4,3);
             zeros(1,4)  zeros(1,4)  -X21          zeros(1,2)  zeros(1,3);
             zeros(2,4)  zeros(2,4)  zeros(2,1)    eye(2)      zeros(2,3);
             zeros(3,4)  zeros(3,4)  zeros(3,1)    zeros(3,2)  -gamma_sqr*eye(3)] * ... 
             PiX2 <= 0;
        PiY2' * ...
            [zeros(4)    Ytt2        zeros(4,1)    zeros(4,2)  zeros(4,3);
             Ytt2        zeros(4)    zeros(4,1)    zeros(4,2)  zeros(4,3);
             zeros(1,4)  zeros(1,4)  -Y21          zeros(1,2)  zeros(1,3);
             zeros(2,4)  zeros(2,4)  zeros(2,1)    eye(2)      zeros(2,3);
             zeros(3,4)  zeros(3,4)  zeros(3,1)    zeros(3,2)  -(1/gamma_sqr)*eye(3)] * ... 
             PiY2 >= 0;
cvx_end 

%% Stop if cvx has failed
if cvx_status ~= 'Solved'
    disp('CVX is not able to solve the problem for this value of gamma.');
    return;
end

%% Check the by cvx computed values:
display('The smallest Xtti eigenvalue is:' + " " + num2str(min([min(eig(Xtt1)), min(eig(Xtt2))])))
display('The smallest Ytti eigenvalue is:' + " " + num2str(min([min(eig(Ytt1)), min(eig(Ytt2))])))


%% Create the controllers
% Just as before; the 1st and Nth pool have some specialties because of the
% interconnection properties. 
X10 = 0; Y10 = 0;

% Function to create the dynamic control matrix of pool Si
S1 = createsi(Ass1, Ast1, Att1, Ats1, Bsn1, Bsu1, Btn1, Btu1, Csy1, Csz1, Cty1, Ctz1, Dyn1, Dzn1, Dzu1, X10, X21, Xtt1, Y10, Y21, Ytt1, gamma_sqr, 1, N);

%% The Nth pool:
% The Nth pool has some pool specific properties in the calculation.
X32 = []; Y32 = []; Ass2 = []; 

% Create SN
S2 = createsi(Ass2, Ast2, Att2, Ats2, Bsn2, Bsu2, Btn2, Btu2, Csy2, Csz2, Cty2, Ctz2, Dyn2, Dzn2, Dzu2, X21, X32, Xtt2, Y21, Y32, Ytt2, gamma_sqr, 2, N);
disp('Si-matrices created');

%% Simulink matrices and models
% Matrices used in the simulink file to split the output of the dynamic
% controller into wiK and uiK
C1 = [1 0]; % w_i^K
C2 = [0 1]; % u_i^K

% Create the state-space A,B,C,D matrices from the computed S matrices
[S1A, S1B, S1C, S1D, ss1] = splitup_S_matrix(S1);
[S2A, S2B, S2C, S2D, ss2] = splitup_S_matrix(S2);

% Now construct discrete-time versions of the dynamic control matrices
h = 1; % Sampling time - STILL NEEDS TO BE PROPERLY SELECTED <------------
% I use the 'tustin' method for phase property preservation of the contr.
ss1d = c2d(ss1, h, 'tustin'); 
ss2d = c2d(ss2, h, 'tustin');

% Compute the delay in number of samples for the discrete-time simulation
for i=1:N
    ddelay(i) = round(tau(i)/h); % Discrete delay
end


% Discretize the shaping weights
W1d = c2d(W1,h,'tustin');
W2d = c2d(W2,h,'tustin');

% Define the continuous-time plant models (third order)
P1 = tf([1],[alpha(1)/w_n(1)^2 2*alpha(1)*zeta/w_n(1) alpha(1) 0 ]);
P2 = tf([1],[alpha(2)/w_n(2)^2 2*alpha(2)*zeta/w_n(2) alpha(2) 0 ]);

% Discretize the plant models
P1d = c2d(P1, h, 'zoh');
P2d = c2d(P2, h, 'zoh');


%% Save the workspace for use in the simulation
% This way, we don't use this script as initialization script, but instead
% we can use a script that simply loads the set of precomputed variables.
save('distributed_workspace.mat'); 
disp(strcat('Problem solved for \gamma = ', num2str(sqrt(gamma_sqr))));
disp('Workspace saved to file for use in simulation.');
































