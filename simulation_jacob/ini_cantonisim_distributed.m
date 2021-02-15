% Author: Jacob Lont
% Date: 05-06-2020
% Initialization script for the simulink file cantonisim_distributed.slx

% Script name: ini_cantonisim_distributed.m

%% Load the precomputed variables
% And throw a custom error message when the set of variables is not found
% and the script that should generate it does not do its job.

try
  load('distributed_workspace.mat');
catch
    cantoni_LMI.m % Run the script that generates the set and try again
    try
        load('distributed_workspace.mat');
    catch
        disp('The set of variables we tried to load is not present. Probably something is wrong with "cantoni_LMI.m"')
    end
end

