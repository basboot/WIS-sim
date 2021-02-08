%% Load data sets for idfentification
% pool1_set1 = wis_data("20210126_step_gate1_2_s25_no_intake.csv", wis);
% pool1_set2 = wis_data("20210126_step_gate1_2_s100_no_intake.csv", wis);
% pool1_set3 = wis_data("20210126_step_gate1_2_s255_no_intake.csv", wis);
% 
% pool2_set1 = wis_data("20210126_step_gate2_3_s25_no_intake.csv", wis);
% pool2_set2 = wis_data("20210126_step_gate2_3_s100_no_intake.csv", wis);
% pool2_set3 = wis_data("20210126_step_gate2_3_s255_no_intake.csv", wis);
% 
% pool3_set1 = wis_data("20210126_step_gate3_4_s25_no_intake.csv", wis);
% pool3_set2 = wis_data("20210126_step_gate3_4_s100_no_intake.csv", wis);
% pool3_set3 = wis_data("20210126_step_gate3_4_s255_no_intake.csv", wis);


%csv_file, wis, pool, type, description, dt

%% Load data sets for identification
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
pool_data(12) = wis_data("20210202_step_gate3_4_s100_no_intake2.csv", wis, 3, "validation", "step 100 pool 3", 1/128);
pool_data(13) = wis_data("20210202_step_gate3_4_s255_no_intake2.csv", wis, 3, "validation", "step 255 pool 3", 1/128);
pool_data(14) = wis_data("20210202_step_gate3_4_s50_no_intake2.csv", wis, 3, "validation", "step 50 pool 3", 1/128);

pool_data(15) = wis_data("20210202_step_gate3_4_s50_no_intake3.csv", wis, 3, "nvalidation", "step 50 pool 3 b", 1/128);

% TODO: create real validation data for pool 1 and 2
pool_data(16) = wis_data("20210202_step_gate1_2_s25_no_intake.csv", wis, 1, "validation", "step 25 pool 1", 1/128);
pool_data(17) = wis_data("20210202_step_gate1_2_s100_no_intake.csv", wis, 1, "validation", "step 100 pool 1", 1/128);
pool_data(18) = wis_data("20210202_step_gate1_2_s255_no_intake.csv", wis, 1, "validation", "step 255 pool 1", 1/128);

pool_data(19) = wis_data("20210202_step_gate2_3_s25_no_intake.csv", wis, 2, "validation", "step 25 pool 2", 1/128);
pool_data(20) = wis_data("20210202_step_gate2_3_s100_no_intake.csv", wis, 2, "validation", "step 100 pool 2", 1/128);
pool_data(21) = wis_data("20210202_step_gate2_3_s255_no_intake.csv", wis, 2, "validation", "step 255 pool 2", 1/128);



