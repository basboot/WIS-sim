%% Define the matrices
% Pool 1
Att1 = 1;
Ats1 = 1;
Btn1 = 1;
Ast1 = zeros(0,1);
Ass1 = zeros(0,1);
Bsn1 = zeros(0,1);
Ctz1 = 1;
Csz1 = 1;
Dzn1 = 0;
Btu1 = 1;
Bsu1 = zeros(0,1);
Dzu1 = 0;

CCD = [Ctz1, Csz1, Dzn1];
BBD = [Btu1', Bsu1', Dzu1'];

N1X = null(CCD);
N1Y = null(BBD);

PiX1 = [1, 0, 0;
        Att1, Ats1, Btn1;
        Ast1, Ass1, Bsn1;
        0, 1, 0;
        Ctz1, Csz1, Dzn1;
        0, 0, 1] * N1X;

PiY1 = [Att1', Ast1', Ctz1';
        -1,   zeros(1,0),    0;
        %0,  -zeros(1,0),  0;       % Only comment for pool 1
        Ats1',  Ass1', Csz1';
        0,    zeros(1,0),    -1;
        Btn1', Bsn1', Dzn1'] * N1Y;
    
% Pool 2

% Pool 3

% Pool 4

% Pool 5


%% Define the optimization problem for the string of pools

cvx_begin sdp  % semi-definite programming
    variable Xtt1(1,1) semidefinite;
    variable Ytt1(1,1) semidefinite;
    %variable X10(nz,nz) symmetric;
    %variable Y10(nz,nz) symmetric;
    variable X21(1,1) symmetric;
    variable Y21(1,1) symmetric;
    variable gamma_sqr nonnegative
    minimize gamma_sqr
    subject to
        Y21 <= 0;
        X21 <= 0;
        [Xtt1, 1; 1, Ytt1] >= 0;
        % [X10 ... ] <= 0
        PiX1' * ...
            [0 Xtt1 0 0 0;
             Xtt1 0 0 0 0;
             0 0 X21 0 0;
             0 0 0 1 0;
             0 0 0 0 -gamma_sqr*1] * PiX1 <= 0;
        PiY1' * ...
            [0 Ytt1 0 0 0;
             Ytt1 0 0 0 0;
             0 0 Y21 0 0;
             0 0 0 1 0;
             0 0 0 0 -gamma_sqr*1] * PiY1 >= 0;
cvx_end 
    
    

