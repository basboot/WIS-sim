%% load_pool_data20210302.m

% new data with more delay

% Load data sets for identification

% All data sets labeled with type "experiment" will be merged for
% identification.
% All data sets labeled "validation" will be used for validation of the
% identification.
% Other types will be ignored.

% 20210302, 128SPS, with extra gates in the pools for delay
PoolData(1) = createWisData("20210302_step_gate1_2_s25_slow_no_intake.csv", Wis, 1, "experiment", "step 25 pool 1", 1/128);
PoolData(2) = createWisData("20210302_step_gate1_2_s100_slow_no_intake.csv", Wis, 1, "validation", "step 100 pool 1", 1/128);
PoolData(3) = createWisData("20210302_step_gate1_2_s255_slow_no_intake.csv", Wis, 1, "experiment", "step 255 pool 1", 1/128);

PoolData(4) = createWisData("20210302_step_gate2_3_s25_slow_no_intake.csv", Wis, 2, "experiment", "step 25 pool 2", 1/128);
PoolData(5) = createWisData("20210302_step_gate2_3_s100_slow_no_intake.csv", Wis, 2, "experiment", "step 100 pool 2", 1/128);
PoolData(6) = createWisData("20210302_step_gate2_3_s255_slow_no_intake.csv", Wis, 2, "validation", "step 255 pool 2", 1/128);

PoolData(7) = createWisData("20210302_step_gate3_4_s25_slow_no_intake.csv", Wis, 3, "experiment", "step 25 pool 3", 1/128);
PoolData(8) = createWisData("20210302_step_gate3_4_s100_slow_no_intake.csv", Wis, 3, "validation", "step 100 pool 3", 1/128);
PoolData(9) = createWisData("20210302_step_gate3_4_s255_slow_no_intake.csv", Wis, 3, "experiment", "step 255 pool 3", 1/128);

