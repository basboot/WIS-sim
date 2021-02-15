%% Construction of a set of distributed controllers for a string of pools
% Source: Distributed controller design for open water channels (2008)
% Yuping Li & Michael Cantoni

% Author: Jacob Lont with help from Gabriel A. Gleizer
% Date: 12-11-2019
% Last modified: 13-11-2019

%% Parameters
clear all; % for debug purposes, can be removed later
% Pool model parameters
alpha = [6492, 2478, 6084, 5658, 7650]; % m^2
tau = [4, 2, 4, 4, 6]; % minutes
phi = [0.48, 1.05, 0.48, 0.48, 0.42 ]; % rad/min

% Loop shaping weights parameters
kappa = [1.69, 6.47, 2.37, 2.21, 1.68];
phi = [113.64, 37.17, 86.96, 96.15, 113.64];
rho = [9.97, 3.26, 7.60, 8.47, 9.97];
eta = [130, 223, 183, 170, 153];

% %% Define the matrices

for i=1:5
    % Create the pool specific matrices for pool i
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

    CCD = [Ctzi, Cszi, Dzni];
    BBD = [Btui', Bsui', Dzui'];

    NiX = null(CCD);
    NiY = null(BBD);

    PiXi = [eye(4),      zeros(4,1),   zeros(4,3);
            Atti,        Atsi,         Btni;
            Asti,        Assi,         Bsni;
            zeros(1,4),  1,            zeros(1,3);
            Ctzi,        Cszi,         Dzni;
            zeros(3,4),  zeros(3,1),   eye(3)]       * NiX;

    PiYi = [Atti',       Asti',        Ctzi';
            -eye(4),     zeros(4,1),   zeros(4,2);
            zeros(1,4),  -eye(1),      zeros(1,2);       % Only comment line for pool 1
            Atsi',       Assi',        Cszi';
            zeros(2,4),  zeros(2,1),   -eye(2);
            Btni',       Bsni',        Dzni']       * NiY;

    switch i
        case 1          % Pool 1 specific matrices
            N1X = NiX;
            % Modify N1Y to comply with the empty (5th) column of PiY1
            N1Y = [NiY(1:4,:);  NiY(6:7,:)]; % Only for pool 1

            Att1 = Atti; Ats1 = Atsi; 
            Btn1 = Btni; Btu1 = Btui;

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

            Att2 = Atti; Ats2 = Atsi; Ast2 = Asti; 
            Btn2 = Btni; Btu2 = Btui;

            Ctz2 = Ctzi; Csz2 = Cszi; Cty2 = Ctyi; Csy2 = Csyi;
            Dzn2 = Dzni; Dzu2 = Dzui; Dyn2 = Dyni; Dyu2 = Dyui;
            
            PiX2 = PiXi;
            PiY2 = PiYi;
            
            
        case 3          % Pool 3 specific matrices
            N3X = NiX;
            N3Y = NiY;

            Att3 = Atti; Ats3 = Atsi; Ast3 = Asti; 
            Btn3 = Btni; Btu3 = Btui;

            Ctz3 = Ctzi; Csz3 = Cszi; Cty3 = Ctyi; Csy3 = Csyi;
            Dzn3 = Dzni; Dzu3 = Dzui; Dyn3 = Dyni; Dyu3 = Dyui;
            
            PiX3 = PiXi;
            PiY3 = PiYi;

        case 4          % Pool 4 specific matrices
            N4X = NiX;
            N4Y = NiY;

            Att4 = Atti; Ats4 = Atsi; Ast4 = Asti; 
            Btn4 = Btni; Btu4 = Btui;

            Ctz4 = Ctzi; Csz4 = Cszi; Cty4 = Ctyi; Csy4 = Csyi;
            Dzn4 = Dzni; Dzu4 = Dzui; Dyn4 = Dyni; Dyu4 = Dyui;
            
            PiX4 = PiXi;
            PiY4 = PiYi;

        case 5          % Pool 5 specific matrices
%             N5X = NiX;
            N5Y = NiY;
            % Modify N5X to comply with the empty (5th) column of PiX5
            % Remove row 5 of N5X
            N5X = [NiX(1:4,:);  NiX(6:8,:)]; % Only for pool 1

            Att5 = Atti; Ast5 = Asti;
            Btn5 = Btni; Btu5 = Btui; Bsn5 = Bsni;

            Ats5 = zeros(4,0); % Only for pool N
            Ass5 = zeros(1,0); % Only for pool N
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
    variable gamma_sqr nonnegative
    %
    minimize gamma_sqr
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
            [zeros(4)    Ytt1       zeros(4,1)   zeros(4,2)  zeros(4,3);
             Ytt1        zeros(4)   zeros(4,1)   zeros(4,2)  zeros(4,3);
             zeros(1,4)  zeros(1,4) Y21          zeros(1,2)  zeros(1,3);
             zeros(2,4)  zeros(2,4)  zeros(2,1)  eye(2)      zeros(2,3);
             zeros(3,4)  zeros(3,4)  zeros(3,1)   zeros(3,2) -(1/gamma_sqr)*eye(3)] * ...
             PiY1 >= 0;
         % 2
         PiX2' * ...
            [0     Xtt2 0    0   0      0;
             Xtt2  0    0    0   0      0;
             0     0    -X21 0   0      0;
             0     0    0    X32 0      0;
             0     0    0    0   eye(2) 0;
             0     0    0    0   0      -gamma_sqr*eye(3)] * PiX2 <= 0;
        PiY2' * ...
            [0    Ytt2 0    0   0      0;
             Ytt2 0    0    0   0      0;
             0    0    -Y21 0   0      0;
             0    0    0    Y32 0      0;
             0    0    0    0   eye(2) 0;
             0    0    0    0   0      -(1/gamma_sqr)*eye(3)] * PiY2 >= 0;
         % 3
         PiX3' * ...
            [0     Xtt3 0    0   0      0;
             Xtt3  0    0    0   0      0;
             0     0    -X32 0   0      0;
             0     0    0    X43 0      0;
             0     0    0    0   eye(2) 0;
             0     0    0    0   0      -gamma_sqr*eye(3)] * PiX3 <= 0;
        PiY3' * ...
            [0    Ytt3 0    0   0      0;
             Ytt3 0    0    0   0      0;
             0    0    -Y32 0   0      0;
             0    0    0    Y43 0      0;
             0    0    0    0   eye(2) 0;
             0    0    0    0   0      -(1/gamma_sqr)*eye(3)] * PiY3 >= 0;
         % 4
         PiX4' * ...
            [0     Xtt4 0    0   0      0;
             Xtt4  0    0    0   0      0;
             0     0    -X43 0   0      0;
             0     0    0    X54 0      0;
             0     0    0    0   eye(2) 0;
             0     0    0    0   0      -gamma_sqr*eye(3)] * PiX4 <= 0;
        PiY4' * ...
            [0    Ytt4 0    0   0      0;
             Ytt4 0    0    0   0      0;
             0    0    -Y43 0   0      0;
             0    0    0    Y54 0      0;
             0    0    0    0   eye(2) 0;
             0    0    0    0   0      -(1/gamma_sqr)*eye(3)] * PiY4 >= 0;
         % 5 (different, Nth pool)
         PiX5' * ...
            [0     Xtt5 0     0      0;
             Xtt5  0    0     0      0;
             0     0    -X54  0      0;
             0     0    0     eye(2) 0;
             0     0    0     0      -gamma_sqr*eye(3)] * PiX5 <= 0;
        PiY5' * ...
            [0    Ytt5 0    0      0;
             Ytt5 0    0    0      0;
             0    0    -Y54 0      0;
             0    0    0    eye(2) 0;
             0    0    0    0      -(1/gamma_sqr)*eye(3)] * PiY5 >= 0;
cvx_end 

%% Check the by cvx computed values:
min(re(eig(Xtt1)),re(eig(Xtt2)),re(eig(Xtt3)),re(eig(Xtt4)), re(eig(Xtt5)))
    

