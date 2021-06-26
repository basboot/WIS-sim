%% WIS (combined plant + local control, cont)

% states of the plant 
%   y_1, water level pool 1
%   delta_1, pade approximation of the delay of pool 1
%   u_1, inflow of pool 1
%   omega_1, loop shaping pole of pool 1
%   y_2,
%   delta_2,
%   u_2,
%   omega_2,
%   y_3,
%   delta_3,
%   u_3,
%   omega_3

Ap = [0 5.39665407447383 -5.39665407447383 0 0 0 -5.39665407447383 0 0 0 0 0 0;
      0 -92.3076923076923 184.615384615385 0 0 0 0 0 0 0 0 0 0;
      0 0 0 1 0 0 0 0 0 0 0 0 0;
      0 0 0 -10 0 0 0 0 0 0 0 0 0;
      0 0 0 0 0 8.424599831508 -8.424599831508 0 0 0 -8.424599831508 0 0;
      0 0 0 0 0 -171.428571428571 342.857142857143 0 0 0 0 0 0;
      0 0 0 0 0 0 0 1 0 0 0 0 0;
      0 0 0 0 0 0 0 -10 0 0 0 0 0;
      0 0 0 0 0 0 0 0 0 4.38788942518649 -4.38788942518649 0 0;
      0 0 0 0 0 0 0 0 0 -80 160 0 0;
      0 0 0 0 0 0 0 0 0 0 0 1 0;
      0 0 0 0 0 0 0 0 0 0 0 -10 0;
      0 0 0 0 0 0 0 0 0 0 0 0 0];  % disturbance

Bp = [0 0 0;0 0 0;30 0 0;-297 0 0;0 0 0;0 0 0;0 50 0;0 -495 0;0 0 0;0 0 0;0 0 30;0 0 -297;
     0 0 0];

Cp = [1 0 0 0 0 0 0 0 0 0 0 0 0;0 0 0 0 1 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 0 1 0 0 0 0];
Cp = [1 0 0 0 0 0 0 0 0 0 0 0 0;
      0 0 0 0 1 0 0 0 0 0 0 0 0;
      0 0 0 0 0 0 0 0 1 0 0 0 0;
      0 0 1 0 0 0 0 0 0 0 0 0 0;
      0 0 0 0 0 0 1 0 0 0 0 0 0;
      0 0 0 0 0 0 0 0 0 0 1 0 0];
      
Dp = [0 0 0;0 0 0;0 0 0];
Dp = 0*Cp*Bp;

C_eye = eye(12);
D_zeros = zeros(12,3);

% noise on pool 3 outflow which is not in the plant 
% results in a dropping water level
Wis.area3 = 0.2279; %m2
% disturbance 0.015 m^3/min
dist = 1 / Wis.area3; % m/min 
E = dist * [0; 0; 0; 0; 0; 0; 0; 0; 1; 0; 0; 0; 0]; 
% This will be an extra, small fluctuation about the steady state
% disturbance.
Ap(9,end) = dist;  % Now in the A matrix.

% WIS controller (P control only)

Ac = zeros(0, 0);
Bc = zeros(0, 3);
Cc = zeros(3, 0);

Dc = [-0.221876653844041 0 0;0 -0.311800712632069 0;0 0 -0.569668315662366];
Bc = zeros(0, 6);
Dc = [Dc, zeros(3)];


% Put time scale to seconds (for conditioning)
Ap = Ap/60;
Bp = Bp/60;
E = E/60;

% Epoch length
h = 1;  % seconds now!


%% Dimensions

np = size(Ap,1);  % states of the plant
nc = size(Ac,1);  % states of the controller
pp = size(Cp,1);  % measured plant outputs
mp = size(Bp,2);  % number of control inputs
nw = size(E,2);   % number of disturbances
ppt = 3;  % For WIS: number of outputs to use in triggering


%% Assumption 4: Bound on disturbance
W_MAG = 0.015;
W_MAG = W_MAG*1000;  % Units to mm and the like
W_MAG = W_MAG/100;  % 1% fluctuation about real deal.

% Bound on noise
if exist('V_GLOBAL','var')
    V_EACH_ELEMENT = V_GLOBAL;
    disp('Using externally set noise value (V_GLOBAL)');
else
    V_EACH_ELEMENT = 0.1; % in mm
end

YFACTOR = 1.1;  % We never know the noise levels that precisely.
V = V_EACH_ELEMENT^2*YFACTOR^2*eye(pp)*pp;