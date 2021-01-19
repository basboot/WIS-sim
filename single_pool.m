% simple first order simulation of a pool
clear all;

s = tf('s');

% pool 0
initial_height0 = 30; % (cm)

% pool 1
tau1 = 0.5;      % delay (s)
alpha1 = 0.1853; % pool area (m^2)


sim('sim_first_order');

% ud = downsample(u, 10);
% yd = downsample(y, 10);
