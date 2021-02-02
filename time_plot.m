%% show graph to check calibration

plot_data = ...
    wis_data("20210202_step_gate3_4_s50_no_intake3.csv", wis);

figure(3);
title("Water heights");
plot(plot_data.timing/1000, plot_data.water_heights);
legend("s1", "s2", "s3", "s4", "s5", "s6", "s7")
                  