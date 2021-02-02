%% show graph to check calibration

%% Load experiment data
step255_2_3 = ...
    wis_data("20210126_step_gate3_4_s255_no_intake.csv", wis);

figure(1);
title("Water heights");
plot(step255_2_3.timing/1000, step255_2_3.water_heights);
legend("s1", "s2", "s3", "s4", "s5", "s6", "s7")
xlabel("time [s]");
ylabel("water height [m]");

figure(2);
[input1, ~] = get_flows_for_pool(step255_2_3.flow_in, step255_2_3.flow_out, 3, true);

output1 = step255_2_3.water_heights(:, 7);
dt1 = step255_2_3.dt;

%% Load validation data

step100_2_3 = ...
    wis_data("20210126_step_gate3_4_s100_no_intake.csv", wis);

figure(3)
[input2, ~] = get_flows_for_pool(step100_2_3.flow_in, step100_2_3.flow_out, 3, true);

output2 = step100_2_3.water_heights(:, 7);
dt2 = step100_2_3.dt;


%% shift data (for delay) and select most reliable part

% Estimating delay not working correctly on this data
% Delays read from step graphs (response time at sensor on other side after 
% opening the gate) => rough estimates!

% pool1 1.2s (20 samples)
% pool2 0.4s (7 samples)
% pool3 1.4s (22 samples)
% 
% [input1, output1] = shift_delay(input1(145:350), output1(145:350), 0);
% [input2, output2] = shift_delay(input2(100:590), output2(100:590), 0);

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
nz = [2];
ioDelay = [22];

mtf = tfest(ze,np,nz,ioDelay,Opt);

%%
figure(8);
compare(zv,mtf)

