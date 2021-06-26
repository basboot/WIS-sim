%% Simulation data
x0 = 0*ones(np,1);  % Assume it's in steady state for the first moment
x0(end) = 100*W_MAG;  % This is the real disturbance
xc0 = zeros(nc,1);
y0 = Cp*x0(1:np);

TEND = 30*60;
kfinal = 30;  % \bar{\kappa} in the paper, heartbeat = 1/2 min

% Disturbance signal
omega = @(t) W_MAG*((t >= 0) & (t <= TEND/2));  % Like in NecSys
odeplant = @(t,xp,u) Ap*xp + Bp*u + E*omega(t);

% For reproducibility, pre-compute the noises
rng(1907);
noises = 2*V_EACH_ELEMENT*(rand(pp,TEND/h + kfinal + 1) - 0.5);

%% Initialize all states need for sleep controller
xp = x0;
y = y0 + noises(:,1);
yhat = y0;
triggered = false;

k = 0;
xc = xc0;
xptilde = zeros(np,1);
X = zeros(np, np);
initialized = false;
psibar = zeros(pp*(kbar+1),1);

sleep = 0;
