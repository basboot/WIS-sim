%% show graph to check calibration

%% Load experiment data
[water_heights, timing, delta_volume, flow_in, flow_out, dt] = ...
    wis_data("20210126_step_gate3_4_s255_no_intake.csv", wis);

figure(1);
title("Water heights");
plot(timing/1000, water_heights);
legend("s1", "s2", "s3", "s4", "s5", "s6", "s7")

figure(2);
[input1, ~] = get_flows_for_pool(flow_in, flow_out, 3, true);

output1 = water_heights(:, 7);
dt1 = dt;

%% Load validation data

[water_heights, timing, delta_volume, flow_in, flow_out, dt] = ...
    wis_data("20210126_step_gate3_4_s100_no_intake.csv", wis);

figure(3)
[input2, ~] = get_flows_for_pool(flow_in, flow_out, 3, true);

output2 = water_heights(:, 7);
dt2 = dt;


%% shift data (for delay) and select most reliable part

% Estimating delay not working correctly on this data
% Delays read from step graphs (response time at sensor on other side after 
% opening the gate) => rough estimates!

% pool1 1.2s (20 samples)
% pool2 0.4s (7 samples)
% pool3 1.4s (22 samples)

[input1, output1] = shift_delay(input1(145:350), output1(145:350), 22);
[input2, output2] = shift_delay(input2(100:590), output2(100:590), 22);

output1 = output1 - output1(1);
output2 = output2 - output2(1);

%% experiment
ze = iddata(output1,input1*1000,dt1); 
figure(4)
plot(ze)

%% validation

zv = iddata(output2,input2*1000,dt2); 
figure(5)
plot(zv)

% %% bode experiment
% 
% Ge = spa(ze);
% figure(5)
% bode(Ge)

%%
% model estimation, experiment
Mimp = impulseest(ze,60);  
figure(6)
% step response
step(Mimp)

%%

% model estimation, validation
Mimp2 = impulseest(zv,60); 
figure(7)
% step response
step(Mimp2)

%% estimate delays
disp('delay')
delayest(ze)
delayest(zv)

%%

Opt = tfestOptions('Display','on');
np = [3]; % 3rd order
nz = [2]
ioDelay = [0];

mtf = tfest(zv,np,nz,ioDelay,Opt);

%%
figure(8);
compare(ze,mtf)

