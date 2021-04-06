%% main.m 

% Runs all scripts needed to identify the lab setup from data and saves 
% the result as a Matlab workspace file.

clear all;

%% Add common functions to path
addpath ../functions/

%% set global properties
wis_properties;

%% perform calibration
calibration;

%% visual validation of calibration
time_plot;

%% load data sets
load_pool_data;
%load_pool_data20210302;


%% identify flow over the gates
gate_identification;

%% identify transfer function for the pools
identification;

%% save results for later use
save("identification.mat", 'Wis', 'PoolModel')

