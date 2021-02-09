%% Load data sets for idfentification

%% Load data sets for identification

% All data sets labeled with type "experiment" will be merged for
% identification.
% All data sets labeled "validation" will be used for validation of the
% identification.
% Other types will be ignored.

% 20210202, 128SPS
pool_data(1) = wis_data("20210202_step_gate1_2_s25_no_intake.csv", wis, 1, "experiment", "step 25 pool 1", 1/128);
pool_data(2) = wis_data("20210202_step_gate1_2_s100_no_intake.csv", wis, 1, "experiment", "step 100 pool 1", 1/128);
pool_data(3) = wis_data("20210202_step_gate1_2_s255_no_intake.csv", wis, 1, "experiment", "step 255 pool 1", 1/128);

pool_data(4) = wis_data("20210202_step_gate2_3_s25_no_intake.csv", wis, 2, "experiment", "step 25 pool 2", 1/128);
pool_data(5) = wis_data("20210202_step_gate2_3_s100_no_intake.csv", wis, 2, "experiment", "step 100 pool 2", 1/128);
pool_data(6) = wis_data("20210202_step_gate2_3_s255_no_intake.csv", wis, 2, "experiment", "step 255 pool 2", 1/128);

pool_data(7) = wis_data("20210202_step_gate3_4_s25_no_intake.csv", wis, 3, "nexperiment", "step 25 pool 3", 1/128);
pool_data(8) = wis_data("20210202_step_gate3_4_s100_no_intake.csv", wis, 3, "experiment", "step 100 pool 3", 1/128);
pool_data(9) = wis_data("20210202_step_gate3_4_s255_no_intake.csv", wis, 3, "experiment", "step 255 pool 3", 1/128);
pool_data(10) = wis_data("20210202_step_gate3_4_s50_no_intake.csv", wis, 3, "experiment", "step 50 pool 3", 1/128);

pool_data(11) = wis_data("20210202_step_gate3_4_s25_no_intake2.csv", wis, 3, "nvalidation", "step 25 pool 3", 1/128);
pool_data(12) = wis_data("20210202_step_gate3_4_s100_no_intake2.csv", wis, 3, "nvalidation", "step 100 pool 3", 1/128);
pool_data(13) = wis_data("20210202_step_gate3_4_s255_no_intake2.csv", wis, 3, "nvalidation", "step 255 pool 3", 1/128);
pool_data(14) = wis_data("20210202_step_gate3_4_s50_no_intake2.csv", wis, 3, "nvalidation", "step 50 pool 3", 1/128);

pool_data(15) = wis_data("20210202_step_gate3_4_s50_no_intake3.csv", wis, 3, "nvalidation", "step 50 pool 3 b", 1/128);

% TODO: Create validation at 128SPS data for pool 1 and 2

pool_data(16) = wis_data("20210202_step_gate1_2_s25_no_intake.csv", wis, 1, "nvalidation", "step 25 pool 1", 1/128);
pool_data(17) = wis_data("20210202_step_gate1_2_s100_no_intake.csv", wis, 1, "nvalidation", "step 100 pool 1", 1/128);
pool_data(18) = wis_data("20210202_step_gate1_2_s255_no_intake.csv", wis, 1, "nvalidation", "step 255 pool 1", 1/128);

pool_data(19) = wis_data("20210202_step_gate2_3_s25_no_intake.csv", wis, 2, "nvalidation", "step 25 pool 2", 1/128);
pool_data(20) = wis_data("20210202_step_gate2_3_s100_no_intake.csv", wis, 2, "nvalidation", "step 100 pool 2", 1/128);
pool_data(21) = wis_data("20210202_step_gate2_3_s255_no_intake.csv", wis, 2, "nvalidation", "step 255 pool 2", 1/128);

% 20210126, 16SPS

pool_data(22) = wis_data("20210126_step_gate1_2_s25_no_intake.csv", wis, 1, "slow", "step 25 pool 1 1/16", 1/16);
pool_data(23) = wis_data("20210126_step_gate1_2_s100_no_intake.csv", wis, 1, "slow", "step 100 pool 1 1/16", 1/16);
pool_data(24) = wis_data("20210126_step_gate1_2_s255_no_intake.csv", wis, 1, "validation", "step 255 pool 1 1/16", 1/16);

pool_data(25) = wis_data("20210126_step_gate2_3_s25_no_intake.csv", wis, 2, "slow", "step 25 pool 1 1/16", 1/16);
pool_data(26) = wis_data("20210126_step_gate2_3_s100_no_intake.csv", wis, 2, "slow", "step 100 pool 1 1/16", 1/16);
pool_data(27) = wis_data("20210126_step_gate2_3_s255_no_intake.csv", wis, 2, "validation", "step 255 pool 1 1/16", 1/16);

pool_data(28) = wis_data("20210126_step_gate3_4_s25_no_intake.csv", wis, 3, "slow", "step 25 pool 1 1/16", 1/16);
pool_data(29) = wis_data("20210126_step_gate3_4_s100_no_intake.csv", wis, 3, "slow", "step 100 pool 1 1/16", 1/16);
pool_data(30) = wis_data("20210126_step_gate3_4_s255_no_intake.csv", wis, 3, "validation", "step 255 pool 1 1/16", 1/16);

