%% time_plot.m

% Show a graph to visualy check the calibration

PlotData = ...
    createWisData("20210202_step_gate2_3_s255_no_intake.csv", ...
    Wis);

figure(3);
title("Water levels");
plot(PlotData.timing/1000, PlotData.water_levels);
legend("s1", "s2", "s3", "s4", "s5", "s6", "s7")
                  