%% show graph to check calibration

[water_heights, timing, delta_volume, flow_in, flow_out, dt] = ...
    wis_data("20210126_step_gate3_4_s25_no_intake.csv", wis);

figure(3);
title("Water heights");
plot(timing/1000, water_heights);
legend("s1", "s2", "s3", "s4", "s5", "s6", "s7")

% data for identification
input1 = flow_out(:,3);
output1 = water_heights(:, 7);
dt1 = dt;
% 
% pool = iddata(output1,input1,dt1);
% 
% clf
% mi = impulseest(pool); % non-parametric (FIR) model
% showConfidence(impulseplot(mi),3); %impulse response with 3 standard
%                                    %deviations confidence region
%                                    
%                                    
%                                    