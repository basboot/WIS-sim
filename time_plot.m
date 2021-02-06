%% show graph to check calibration

plot_data = ...
    wis_data("20210202_step_gate2_3_s255_no_intake.csv", wis);

figure(3);
title("Water levels");
plot(plot_data.timing/1000, plot_data.water_levels);
legend("s1", "s2", "s3", "s4", "s5", "s6", "s7")
                  