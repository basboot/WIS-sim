%% Construction of a set of distributed controllers for a string of pools
% Source: Distributed controller design for open water channels (2008) [1]
% Yuping Li & Michael Cantoni

% Author: Jacob Lont with help from Gabriel A. Gleizer
% Date: 12-11-2019
% Last modified: 18-06-2020
% Status: Controllers work, but the control may be improved based on 
% comparison with a PI-controller from literature.


%% Parameters
clear all; % for debug purposes, can be removed later

%% Load identification results
try
    load('../identification/identification.mat');
catch
    assert(false, "File 'identification.mat' does not exist. Run identification first.");
end

% Pool model parameters
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


% 3rd order model is used for simulation
% w_n = [0.65, 0.86, 0.448, 0.448, 0.448]; Wrong: computed wrongly
% We can just compute w_n now that we have a sound expression:
zeta = 0.0151;  % no unit % From Gabriel's ini script
w_n = zeros(1,length(phi_wave));
for i=1:length(phi_wave)
    w_n(i) = phi_wave(i)/sqrt(1-zeta^2); % Literature survey sec. 3-3-1.
%     display(w_n(i));
end



% Loop shaping weights parameters from [1]
kappa = [1.69, 6.47, 2.37, 2.21, 1.68];
phi = [113.64, 37.17, 86.96, 96.15, 113.64];
rho = [9.97, 3.26, 7.60, 8.47, 9.97];
eta = [130, 223, 183, 170, 153];

% %% Define the matrices

for i=1:5
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
%     CCD = [Ctzi, Cszi, Dzni]; % <-- 23-07-2020: this is not what the paper has
    CCD2 = [Ctyi, Csyi, Dyni]; % <-- 23-07-2020: This is what the paper has
    BBD = [Btui', Bsui', Dzui'];

%     NiX = null(CCD); % <--- 23-07-2020 Changed this to use CCD2, result
%     is the same matrix NiX.
    NiX = null(CCD2);
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

        case 2          % Pool 2 specific matrices
            N2X = NiX;
            N2Y = NiY;

            Att2 = Atti; Ats2 = Atsi; Ast2 = Asti; Ass2 = Assi;
            Btn2 = Btni; Btu2 = Btui; Bsu2 = Bsui; Bsn2 = Bsni;

            Ctz2 = Ctzi; Csz2 = Cszi; Cty2 = Ctyi; Csy2 = Csyi;
            Dzn2 = Dzni; Dzu2 = Dzui; Dyn2 = Dyni; Dyu2 = Dyui;
            
            PiX2 = PiXi;
            PiY2 = PiYi;
            
            
        case 3          % Pool 3 specific matrices
            N3X = NiX;
            N3Y = NiY;

            Att3 = Atti; Ats3 = Atsi; Ast3 = Asti; Ass3 = Assi;
            Btn3 = Btni; Btu3 = Btui; Bsu3 = Bsui; Bsn3 = Bsni;

            Ctz3 = Ctzi; Csz3 = Cszi; Cty3 = Ctyi; Csy3 = Csyi;
            Dzn3 = Dzni; Dzu3 = Dzui; Dyn3 = Dyni; Dyu3 = Dyui;
            
            PiX3 = PiXi;
            PiY3 = PiYi;

        case 4          % Pool 4 specific matrices
            N4X = NiX;
            N4Y = NiY;

            Att4 = Atti; Ats4 = Atsi; Ast4 = Asti; Ass4 = Assi;
            Btn4 = Btni; Btu4 = Btui; Bsu4 = Bsui; Bsn4 = Bsni;

            Ctz4 = Ctzi; Csz4 = Cszi; Cty4 = Ctyi; Csy4 = Csyi;
            Dzn4 = Dzni; Dzu4 = Dzui; Dyn4 = Dyni; Dyu4 = Dyui;
            
            PiX4 = PiXi;
            PiY4 = PiYi;

        case 5          % Pool 5 specific matrices
%             N5X = NiX;
            N5Y = NiY;
            % Modify N5X to comply with the empty (5th) column of PiX5
            % Remove row 5 of N5X
            N5X = [NiX(1:4,:);  NiX(6:8,:)]; % Only for pool 5

            Att5 = Atti; Ast5 = Asti; %Ass5 = Assi;
            Btn5 = Btni; Btu5 = Btui; Bsn5 = Bsni; Bsu5 = Bsui;

            Ats5 = zeros(4,0); % Only for pool N
            Ass5 = zeros(1,0); % Only for pool N
%             Ass5 = []; % This is what we use later on to Create S5, so
%             maybe we should also implement this here?
            Csz5 = zeros(2,0); % Only for pool N
            
            Ctz5 = Ctzi; Cty5 = Ctyi; Csy5 = Csyi;
            Dzn5 = Dzni; Dzu5 = Dzui; Dyn5 = Dyni; Dyu5 = Dyui;

            % PiX5, PiY5 computation differs from rest of the pools
            PiX5 = [eye(4),     zeros(4,0),  zeros(4,3);
                    Att5,       Ats5,        Btn5;
                    Ast5,       Ass5,        Bsn5;
                    %zeros(1,4), 1,           zeros(1,3); % Only pool 5
                    Ctz5,       Csz5,        Dzn5;
                    zeros(3,4), zeros(3,0),  eye(3)]      * N5X;
            
            PiY5 = [Att5',      Ast5',        Ctz5';
                    -eye(4),    zeros(4,1),   zeros(4,2);
                    zeros(1,4)  -eye(1),      zeros(1,2);
                    Ats5',      Ass5',        Csz5';
                    zeros(2,4), zeros(2,1),   -eye(2);
                    Btn5',      Bsn5',        Dzn5']      * N5Y;
    end

end %for







%% Define the optimization problem for the string of pools

% In case of errors
% Error: Undefined function 'vec' for input arguments of type 'double'.
% Run:
% addpath([cvx_where,'functions',filesep,'vec_']);
% savepath

% gamma_sqr = 8.2944; % <= gamma = 2.88 is optimum from the paper
% Smallest that can be solved: 
% gamma_sqr = 16;
gamma_sqr = 15.4;
disp(strcat('solving for \gamma = ', num2str(sqrt(gamma_sqr))));



cvx_begin sdp  % semi-definite programming
    variable Xtt1(4,4) semidefinite; % implies pos def
    variable Xtt2(4,4) semidefinite;
    variable Xtt3(4,4) semidefinite;
    variable Xtt4(4,4) semidefinite;
    variable Xtt5(4,4) semidefinite;
    variable Ytt1(4,4) semidefinite;
    variable Ytt2(4,4) semidefinite;
    variable Ytt3(4,4) semidefinite;
    variable Ytt4(4,4) semidefinite;
    variable Ytt5(4,4) semidefinite;
    %
    variable X21(1,1) symmetric; % symmetric, but not pos def
    variable X32(1,1) symmetric;
    variable X43(1,1) symmetric;
    variable X54(1,1) symmetric;
    variable Y21(1,1) symmetric;
    variable Y32(1,1) symmetric;
    variable Y43(1,1) symmetric;
    variable Y54(1,1) symmetric;
%     variable gamma_sqr nonnegative
    %
    minimize 0  
%     minimize gamma_sqr
    subject to
        % Did not include: Xitt > 0, Yitt > 0, not sure if this is already
        % implicitly declared or not. Check this later on.
        Y21 <= 0; Y32 <= 0; Y43 <= 0; Y54 <= 0; % < 0
        X21 <= 0; X32 <= 0; X43 <= 0; X54 <= 0; % < 0
        % 
        [Xtt1, eye(4); eye(4), Ytt1] >= 0;
        [Xtt2, eye(4); eye(4), Ytt2] >= 0;
        [Xtt3, eye(4); eye(4), Ytt3] >= 0;
        [Xtt4, eye(4); eye(4), Ytt4] >= 0;
        [Xtt5, eye(4); eye(4), Ytt5] >= 0;
        % 
        [X21, -1; -1, Y21] <= 0;
        [X32, -1; -1, Y32] <= 0;
        [X43, -1; -1, Y43] <= 0;
        [X54, -1; -1, Y54] <= 0;
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
         % 2
         PiX2' * ...
            [zeros(4)    Xtt2        zeros(4,1)   zeros(4,1) zeros(4,2)  zeros(4,3);
             Xtt2        zeros(4)    zeros(4,1)   zeros(4,1) zeros(4,2)  zeros(4,3);
             zeros(1,4)  zeros(1,4)  -X21         0          zeros(1,2)  zeros(1,3);
             zeros(1,4)  zeros(1,4)  0            X32        zeros(1,2)  zeros(1,3);
             zeros(2,4)  zeros(2,4)  zeros(2,1)   zeros(2,1) eye(2)      zeros(2,3);
             zeros(3,4)  zeros(3,4)  zeros(3,1)   zeros(3,1) zeros(3,2)  -gamma_sqr*eye(3)] * ... 
             PiX2 <= 0;
        PiY2' * ...
            [zeros(4)    Ytt2        zeros(4,1)   zeros(4,1) zeros(4,2)  zeros(4,3);
             Ytt2        zeros(4)    zeros(4,1)   zeros(4,1) zeros(4,2)  zeros(4,3);
             zeros(1,4)  zeros(1,4)  -Y21         0          zeros(1,2)  zeros(1,3);
             zeros(1,4)  zeros(1,4)  0            Y32        zeros(1,2)  zeros(1,3);
             zeros(2,4)  zeros(2,4)  zeros(2,1)   zeros(2,1) eye(2)      zeros(2,3);
             zeros(3,4)  zeros(3,4)  zeros(3,1)   zeros(3,1) zeros(3,2)  -(1/gamma_sqr)*eye(3)] * ... 
             PiY2 >= 0;
         % 3
         PiX3' * ...
            [zeros(4)    Xtt3        zeros(4,1)   zeros(4,1) zeros(4,2)  zeros(4,3);
             Xtt3        zeros(4)    zeros(4,1)   zeros(4,1) zeros(4,2)  zeros(4,3);
             zeros(1,4)  zeros(1,4)  -X32         0          zeros(1,2)  zeros(1,3);
             zeros(1,4)  zeros(1,4)  0            X43        zeros(1,2)  zeros(1,3);
             zeros(2,4)  zeros(2,4)  zeros(2,1)   zeros(2,1) eye(2)      zeros(2,3);
             zeros(3,4)  zeros(3,4)  zeros(3,1)   zeros(3,1) zeros(3,2)  -gamma_sqr*eye(3)] * ... 
             PiX3 <= 0;
        PiY3' * ...
            [zeros(4)    Ytt3        zeros(4,1)   zeros(4,1) zeros(4,2)  zeros(4,3);
             Ytt3        zeros(4)    zeros(4,1)   zeros(4,1) zeros(4,2)  zeros(4,3);
             zeros(1,4)  zeros(1,4)  -Y32         0          zeros(1,2)  zeros(1,3);
             zeros(1,4)  zeros(1,4)  0            Y43        zeros(1,2)  zeros(1,3);
             zeros(2,4)  zeros(2,4)  zeros(2,1)   zeros(2,1) eye(2)      zeros(2,3);
             zeros(3,4)  zeros(3,4)  zeros(3,1)   zeros(3,1) zeros(3,2)  -(1/gamma_sqr)*eye(3)] * ... 
             PiY3 >= 0;
         % 4
         PiX4' * ...
            [zeros(4)    Xtt4        zeros(4,1)   zeros(4,1) zeros(4,2)  zeros(4,3);
             Xtt4        zeros(4)    zeros(4,1)   zeros(4,1) zeros(4,2)  zeros(4,3);
             zeros(1,4)  zeros(1,4)  -X43         0          zeros(1,2)  zeros(1,3);
             zeros(1,4)  zeros(1,4)  0            X54        zeros(1,2)  zeros(1,3);
             zeros(2,4)  zeros(2,4)  zeros(2,1)   zeros(2,1) eye(2)      zeros(2,3);
             zeros(3,4)  zeros(3,4)  zeros(3,1)   zeros(3,1) zeros(3,2)  -gamma_sqr*eye(3)] * ... 
             PiX4 <= 0;
        PiY4' * ...
            [zeros(4)    Ytt4        zeros(4,1)   zeros(4,1) zeros(4,2)  zeros(4,3);
             Ytt4        zeros(4)    zeros(4,1)   zeros(4,1) zeros(4,2)  zeros(4,3);
             zeros(1,4)  zeros(1,4)  -Y43         0          zeros(1,2)  zeros(1,3);
             zeros(1,4)  zeros(1,4)  0            Y54        zeros(1,2)  zeros(1,3);
             zeros(2,4)  zeros(2,4)  zeros(2,1)   zeros(2,1) eye(2)      zeros(2,3);
             zeros(3,4)  zeros(3,4)  zeros(3,1)   zeros(3,1) zeros(3,2)  -(1/gamma_sqr)*eye(3)] * ... 
             PiY4 >= 0;
         % 5 (different, Nth pool)
         PiX5' * ...
            [zeros(4)    Xtt5        zeros(4,1)    zeros(4,2)  zeros(4,3);
             Xtt5        zeros(4)    zeros(4,1)    zeros(4,2)  zeros(4,3);
             zeros(1,4)  zeros(1,4)  -X54          zeros(1,2)  zeros(1,3);
             zeros(2,4)  zeros(2,4)  zeros(2,1)    eye(2)      zeros(2,3);
             zeros(3,4)  zeros(3,4)  zeros(3,1)    zeros(3,2)  -gamma_sqr*eye(3)] * ... 
             PiX5 <= 0;
        PiY5' * ...
            [zeros(4)    Ytt5        zeros(4,1)    zeros(4,2)  zeros(4,3);
             Ytt5        zeros(4)    zeros(4,1)    zeros(4,2)  zeros(4,3);
             zeros(1,4)  zeros(1,4)  -Y54          zeros(1,2)  zeros(1,3);
             zeros(2,4)  zeros(2,4)  zeros(2,1)    eye(2)      zeros(2,3);
             zeros(3,4)  zeros(3,4)  zeros(3,1)    zeros(3,2)  -(1/gamma_sqr)*eye(3)] * ... 
             PiY5 >= 0;
cvx_end 

%% Stop if cvx has failed
if strcmp(cvx_status, 'Failed')
    disp('CVX is not able to solve the problem for this value of gamma.');
    return;
end

%% Check the by cvx computed values:
display('The smallest Xtti eigenvalue is:' + " " + num2str(min([min(eig(Xtt1)), min(eig(Xtt2)), min(eig(Xtt3)), min(eig(Xtt4)), min(eig(Xtt5))])))
display('The smallest Ytti eigenvalue is:' + " " + num2str(min([min(eig(Ytt1)), min(eig(Ytt2)), min(eig(Ytt3)), min(eig(Ytt4)), min(eig(Ytt5))])))



%% Create the controllers
% Just as before; the 1st and Nth pool have some specialties because of the
% interconnection properties. 
N = 5; % The number of pools. To show the 'createsi' function which pool is the Nth pool.
% X10 = []; Y10 = []; This does not work out as Xiim1C gets all empty
% entries, which makes it not possible to create Gi as it looks like.
X10 = 0; Y10 = 0;

% clc;
% Function to create the dynamic control matrix of pool 1: S1
S1 = createsi(Ass1, Ast1, Att1, Ats1, Bsn1, Bsu1, Btn1, Btu1, Csy1, Csz1, Cty1, Ctz1, Dyn1, Dzn1, Dzu1, X10, X21, Xtt1, Y10, Y21, Ytt1, gamma_sqr, 1, N);

%% The middle pools. All fine. No exceptions.
S2 = createsi(Ass2, Ast2, Att2, Ats2, Bsn2, Bsu2, Btn2, Btu2, Csy2, Csz2, Cty2, Ctz2, Dyn2, Dzn2, Dzu2, X21, X32, Xtt2, Y21, Y32, Ytt2, gamma_sqr, 2, N);
S3 = createsi(Ass3, Ast3, Att3, Ats3, Bsn3, Bsu3, Btn3, Btu3, Csy3, Csz3, Cty3, Ctz3, Dyn3, Dzn3, Dzu3, X32, X43, Xtt3, Y32, Y43, Ytt3, gamma_sqr, 3, N);
S4 = createsi(Ass4, Ast4, Att4, Ats4, Bsn4, Bsu4, Btn4, Btu4, Csy4, Csz4, Cty4, Ctz4, Dyn4, Dzn4, Dzu4, X43, X54, Xtt4, Y43, Y54, Ytt4, gamma_sqr, 4, N);

%% The Nth pool:
% The Nth pool has some pool specific properties in the calculation.
X65 = []; % Set to zero. Should be a zero dimension entry. It does not exist.
Y65 = [];
Ass5 = [];

% Create S5
S5 = createsi(Ass5, Ast5, Att5, Ats5, Bsn5, Bsu5, Btn5, Btu5, Csy5, Csz5, Cty5, Ctz5, Dyn5, Dzn5, Dzu5, X54, X65, Xtt5, Y54, Y65, Ytt5, gamma_sqr, 5, N);
disp('Si-matrices created');

%% Simulink matrices and models
% Matrices used in the simulink file to split the output of the dynamic
% controller into wiK and uiK
C1 = [1 0]; % w_i^K
C2 = [0 1]; % u_i^K

% Create the state-space A,B,C,D matrices from the computed S matrices
[S1A, S1B, S1C, S1D, ss1] = splitup_S_matrix(S1);
[S2A, S2B, S2C, S2D, ss2] = splitup_S_matrix(S2);
[S3A, S3B, S3C, S3D, ss3] = splitup_S_matrix(S3);
[S4A, S4B, S4C, S4D, ss4] = splitup_S_matrix(S4);
[S5A, S5B, S5C, S5D, ss5] = splitup_S_matrix(S5);

% Now construct discrete-time versions of the dynamic control matrices
h = 1; % Sampling time - STILL NEEDS TO BE PROPERLY SELECTED <------------
% I use the 'tustin' method for phase property preservation of the contr.
[ss1d, ss1d_map] = c2d(ss1, h, 'tustin'); 
[ss2d, ss2d_map] = c2d(ss2, h, 'tustin');
[ss3d, ss3d_map] = c2d(ss3, h, 'tustin');
[ss4d, ss4d_map] = c2d(ss4, h, 'tustin');
[ss5d, ss5d_map] = c2d(ss5, h, 'tustin');

% Compute the delay in number of samples for the discrete-time simulation
for i=1:5
    ddelay(i) = round(tau(i)/h); % Discrete delay
end

% Continuous-time shaping weights
W1 = tf([kappa(1)*phi(1) kappa(1)], [rho(1) 1 0]);
W2 = tf([kappa(2)*phi(2) kappa(2)], [rho(2) 1 0]);
W3 = tf([kappa(3)*phi(3) kappa(3)], [rho(3) 1 0]);
W4 = tf([kappa(4)*phi(4) kappa(4)], [rho(4) 1 0]);
W5 = tf([kappa(5)*phi(5) kappa(5)], [rho(5) 1 0]);

% Discretize the shaping weights
W1d = c2d(W1,h,'tustin');
W2d = c2d(W2,h,'tustin');
W3d = c2d(W3,h,'tustin');
W4d = c2d(W4,h,'tustin');
W5d = c2d(W5,h,'tustin');

% Define the continuous-time plant models (third order)
P1 = tf([1],[alpha(1)/w_n(1)^2 2*alpha(1)*zeta/w_n(1) alpha(1) 0 ]);
P2 = tf([1],[alpha(2)/w_n(2)^2 2*alpha(2)*zeta/w_n(2) alpha(2) 0 ]);
P3 = tf([1],[alpha(3)/w_n(3)^2 2*alpha(3)*zeta/w_n(3) alpha(3) 0 ]);
P4 = tf([1],[alpha(4)/w_n(4)^2 2*alpha(4)*zeta/w_n(4) alpha(4) 0 ]);
P5 = tf([1],[alpha(5)/w_n(5)^2 2*alpha(5)*zeta/w_n(5) alpha(5) 0 ]);

% Discretize the plant models
P1d = c2d(P1, h, 'zoh');
P2d = c2d(P2, h, 'zoh');
P3d = c2d(P3, h, 'zoh');
P4d = c2d(P4, h, 'zoh');
P5d = c2d(P4, h, 'zoh');



%% Save the workspace for use in the simulation (moved to later on)
% This way, we don't use this script as initialization script, but instead
% we can use a script that simply loads the set of precomputed variables.
% save('distributed_workspace.mat'); 
disp(strcat('Problem solved for \gamma = ', num2str(sqrt(gamma_sqr))));
% disp('Workspace saved to file for use in simulation.');


%% Check stability using the LMI script based on Heemels(2013) PETC for ..

gamma = sqrt(gamma_sqr);

% System matrices using the generalized plant defined in eq. 4 in Cantoni
% et al. (2008) - Distributed controller design for open water channels

% Bsu1 = zeros(0,1); % This is different from before but this is right!!
% But, it gives problems, as we now do not have a square A matrix
% Adjust the matrices of the wi row of G1 to be 1 row, 0-valued entries
% Bsu1 = zeros(1,1);
% Bsn1 = zeros(1,3);
% Ast1 = zeros(1,4);
% Ass1 = zeros(1,1);
% Same problems for the number N(=5) case, but now for the Atsi column.
% Ats5 = zeros(4,1);
% Ass5 = zeros(1,1);
% Csz5 = zeros(2,1);
% Csy5 = zeros(1,1);
% --> Chosen to define A in a different way -> A = Att; instead of saying 
% A = [Att, Ats; Ast, Ass]; 

% TODO: Adjust CreateSi to handle this and check!
% Csy5 = zeros(1,0); % This is different from before, but this is right!


% Try to build the interconnection using append
% G = append(G1, G2, G3, G4, G5); % Does give Ap, but the interconnection is not reflected now..

% Now define the combined A,B,C,D plant-matrices
Ap = [Att1, zeros(4,2), Ats1, zeros(4,13);
      zeros(4,4), Att2, zeros(4,2), Ats2, zeros(4,9);
      zeros(4,8), Att3, zeros(4,2), Ats3, zeros(4,5);
      zeros(4,12), Att4, zeros(4,2), Ats4, zeros(4,1);
      zeros(4,16), Att5];
Bp = blkdiag(Btu1, Btu2, Btu3, Btu4, Btu5);

Btnw1 = Btn1(:,2:3); 
Btnw2 = Btn2(:,2:3);
Btnw3 = Btn3(:,2:3);
Btnw4 = Btn4(:,2:3);
Btnw5 = Btn5(:,2:3);
Bw = blkdiag(Btnw1, Btnw2, Btnw3, Btnw4, Btnw5);

Cpi = [1 0 0 0];
Cp = blkdiag(Cpi, Cpi, Cpi, Cpi, Cpi);

% Define the controller by combining the individual controller matrices
% I have the discrete time state-space description of the controller
% Now extend it to have w_i^K as a state.
% Create the submatrices I need, to redefine my states
[Ac1, Bc1, Cc1, Dc1] = create_combined_control_sub_matrices(ss1d, 1);
[Ac2, Bc2, Cc2, Dc2] = create_combined_control_sub_matrices(ss2d, 2);
[Ac3, Bc3, Cc3, Dc3] = create_combined_control_sub_matrices(ss3d, 3);
[Ac4, Bc4, Cc4, Dc4] = create_combined_control_sub_matrices(ss4d, 4);
[Ac5, Bc5, Cc5, Dc5] = create_combined_control_sub_matrices(ss5d, 5);

Ac = [Ac1; Ac2; Ac3; Ac4; Ac5]; % concatenating block rows
Bc = blkdiag(-Bc1, Bc2, Bc3, Bc4, Bc5);
Cc = [Cc1; Cc2; Cc3; Cc4; Cc5]; % concatenating block rows
Dc = blkdiag(-Dc1, Dc2, Dc3, Dc4, Dc5);

nDc = size(Dc,1);
D = [zeros(nDc), zeros(nDc); Dc, zeros(nDc)];
nD = size(D,1); % D should be square (10x10), check eq. 43 of Heemels(2013)
C = blkdiag(Cp,Cc);

%% Create state-space models from the matrices for simulation checks
% Create a simulation to check if the matrix operations are right
% All discrete-time!
Dp = zeros(size(Cp,1),size(Bp,2)); % Because we dont have a Dp
comb_plant_cont = ss(Ap, Bp, Cp, Dp); % combined plant model cont.-time
comb_plant_disc = c2d(comb_plant_cont,h,'zoh'); % combined plant model disc.-time
comb_contr = ss(Ac, Bc, Cc, Dc, h); % combined controller model (discrete-time)

% Simulation works!

%% Check nominal stability using feedback and eigenvalues within unit circle
% They are all inside the unit circle! => Nominal stability
CPcomb = feedback(comb_plant_disc*comb_contr,eye(5),+1); % Combined system, pos. feedback
poles = eig(CPcomb);
fprintf('The eigenvalues of the CL-system are between %f and %f \n',min(real(poles)), max(poles));

if max(abs(poles))>1
    fprintf('NO nominal stability!\n Check the closed-loop poles.\n');
else
    fprintf('Nominal stability: all poles inside the unit circle.\n');
end

% Plot the eigenvalues with the unit circle for visual inspection
% figure(1); clf;
% plot(poles,'*')
% hold on
% ezplot(@(x,y) (x).^2 + (y).^2 -1^2)
% xlim([-1.1,1.1]); ylim([-1.1,1.1]);
% title('Closed loop eigenvalues and the unit circle');
% xlabel('Re'); ylabel('Im');
% hold off

% Compute a reference for rho I put into the LMI:
a = max(abs(poles));
rho_reference = -log(abs(a))/h; % 0.0044

%% Do the Heemels LMI check for GES and L2-gain
if gamma < sqrt(max(eig(D'*D))) %check of Heemels's Theorem V.2, eq 48
    fprintf('gamma is too small according to Heemels for this matrix D.\n');
    % Then M is not invertible see eq. 20 of Heemels 2013
end

% TODO: renamed rho to rho_lmi because it shadows rho used in the
%       simulation => CHECK WHERE RHO IS USED!

% Triggering parameters
rho_lmi = 0.0000000000001; %Tuned manually % rho > 0 , (lower bound on) the decay rate
% sigma = 0.15; % Tuned manually together with rho and the L2gain: 0.15 works nice in simulation
sigma = 0.10;
% sigma = 0.05; % 0.05 Very nice performance and also a great reduction of communication!
% sigma = 0.01; % sigma for LMI: 0.01 looks great, and still less than 25% of samples used
lambda = 20000000; % called gamma in Heemels 2013. L2 upper bound.

% Create Q
% Gamma = eye(10); % In our case of 1 artificial node, based: on 2*y, 2*u
% ng = size(Gamma, 1);
Q11 = (1-sigma)*C'*C;
Q12 = (1-sigma)*C'*D-C';
Q21 = (1-sigma)* D'*C - C;
Q22 = (D-eye(nD))'*(D-eye(nD)) - sigma*D'*D;
Q = [Q11, Q12; Q21, Q22];


%% Save the workspace for use in the simulation
% This way, we don't use this script as initialization script, but instead
% we can use a script that simply loads the set of precomputed variables.
save('distributed_workspace.mat'); 
disp(strcat('Problem solved for \gamma = ', num2str(sqrt(gamma_sqr))));
disp('Workspace saved to file for use in simulation.');



%% Run the LMI solver to check GES and get the L2-gain <= gamma
% Uncomment the next 3 lines to run the Heemels LMI check for output feeedback
% clc;
% fprintf('Calling the LMI-check-function.\n');
% [status, Ph, mu] = checkheemelslmiOF(comb_contr.A, comb_contr.B, comb_contr.C, comb_contr.D, Ap, Bp, Bw, Cp, Q, h, rho_lmi, lambda);
% 





% Without the second LMI added for J0, it solved for 
% sigma = 0.01;  rho = 0.0000000000001; lambda = 2000000000000;
% Also with the second LMI added, but the matrix of the second LMI has a
% complex value for min(eig) and even a negative real eig -5.2003e-05
% this is also achieved for smaller lambda (around 200(0) or so -> check)
% and rho around 0.000001

% Gabriel: Be very careful in the wording about stability. 'numerical
% issues', 'it is likely that we have stability, although because of the
% very numerical instability' we have not been able to get an actual
% positive definite solution.
% Mention the large size of the LMIs, and the ill conditioned matrices:
% very large and very small matrices

% You could try to use cvx_solver sedumi as an alternative -> fails all the
% time

% mu1 and mu2 are independent now. At first we only had 1 mu

% rho = 0.0000001; sigma = 0.01; lambda = 200000; FAILED 
% with: precision high, and the >= 0 conditions

% rho = 0.0000000000001; sigma = 0.01; lambda = 2000000000000; FAILED
% with: precision high, and the >= 0 conditions

% rho = 0.0000000000001; sigma = 0.01; lambda = 2000000000000; FAILED
% with: precision standard, and the >= 0 conditions

% rho = 0.0000000000001; sigma = 0.01; lambda = 2000000000000; FAILED
% with: precision standard, and NOT the >= 0 conditions for Ph (choosing semidefinite),
% but using the >= 0 conditions for mu1 and mu2

% rho = 0.0000000000001; sigma = 0.01; lambda = 2000000000000; FAILED
% with: precision standard, and NOT the >= 0 conditions for Ph (choosing semidefinite),
% and NOT using the >= 0 conditions for mu1 and mu2, but choosing
% nonnegative instead for both mu1 and mu2.

% rho = 0.0000000000001; sigma = 0.01; lambda = 2000000000000; FAILED but
% it was close. It looked like it was converging.
% with: precision standard, USING ONLY mu1 (in both LMIs)
% and NOT the >= 0 conditions for Ph (choosing semidefinite),
% and NOT using the >= 0 conditions for mu1 and mu2, but choosing
% nonnegative instead for both mu1 and mu2.

% Idea: Maybe I should go back to 1 single mu instad of mu1 and mu2,
% however I did get some results at some point using mu1 and mu2 ??
% OUTCOME: It doesnt directly solve my problem.

% rho = 0.0000001; sigma = 0.01; lambda = 2000;  FAILED
% with: precision standard, USING ONLY mu1 (in both LMIs)
% and NOT the >= 0 conditions for Ph (choosing semidefinite),
% and NOT using the >= 0 conditions for mu1 and mu2, but choosing
% nonnegative instead for both mu1 and mu2.

% rho = 0.0000001; sigma = 0.01; lambda = 2000; FAILED
% with: precision HIGH, using both mu1 and mu2
% and NOT the >= 0 conditions for Ph (choosing semidefinite),
% and NOT using the >= 0 conditions for mu1 and mu2, but choosing
% nonnegative instead for both mu1 and mu2.

% rho = 0.0000001; sigma = 0.01; lambda = 20000;   FAILED
% with: precision HIGH, using both mu1 and mu2
% and NOT the >= 0 conditions for Ph (choosing semidefinite),
% and NOT using the >= 0 conditions for mu1 and mu2, but choosing
% nonnegative instead for both mu1 and mu2.

% Idea: Change the performance output matrix Cbar as Gabriel suggested
% -> Changed the matrix Cbar to use the actual outputs instad of y_hat

% rho = 0.0000001; sigma = 0.01; lambda = 200000;   FAILED??
% with: precision HIGH, using both mu1 and mu2
% and NOT the >= 0 conditions for Ph (choosing semidefinite),
% and NOT using the >= 0 conditions for mu1 and mu2, but choosing
% nonnegative instead for both mu1 and mu2.

% rho = 0.0000001; sigma = 0.01; lambda = 200000;   FAILED
% with: precision DEFAULT, using both mu1 and mu2
% and NOT the >= 0 conditions for Ph (choosing semidefinite),
% and NOT using the >= 0 conditions for mu1 and mu2, but choosing
% nonnegative instead for both mu1 and mu2

% rho = 0.0000000000001; sigma = 0.01; lambda = 200000;   FAILED
% with: precision DEFAULT, using both mu1 and mu2
% and NOT the >= 0 conditions for Ph (choosing semidefinite),
% and NOT using the >= 0 conditions for mu1 and mu2, but choosing
% nonnegative instead for both mu1 and mu2

% rho = 0.0000000000001; sigma = 0.001; lambda = 200000;   FAILED
% with: precision DEFAULT, using both mu1 and mu2
% and NOT the >= 0 conditions for Ph (choosing semidefinite),
% and NOT using the >= 0 conditions for mu1 and mu2, but choosing
% nonnegative instead for both mu1 and mu2

% rho = 0.0000000000001; sigma = 0.0001; lambda = 20000000;   FAILED
% with: precision LOW, using both mu1 and mu2
% and NOT the >= 0 conditions for Ph (choosing semidefinite),
% and NOT using the >= 0 conditions for mu1 and mu2, but choosing
% nonnegative instead for both mu1 and mu2












