clear all;

% set global properties
wis_properties;

% perform calibration
calibration;

% visual validation of calibration
time_plot;

% load data sets
load_pool_data;

% identify flow over the gates
gate_identification;

% identify transfer function for the pools
identification;

