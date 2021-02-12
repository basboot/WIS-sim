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

%% identify flow over the gates
gate_identification;

%% identify transfer function for the pools
identification;

%% save results for later use
save("identification.mat", 'wis', 'pool_model')

