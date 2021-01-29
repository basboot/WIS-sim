%% show graph to check calibration

% TODO: Cleanup these tests
% but for an overview:
% inputXa = flow caculated based on pools before
% inputXb = flow caculated based on pools after
% inputX = flow estimated by hand (simple linear)
% outputXa = measured y
% outputX = simple 1st order data




[water_heights, timing, delta_volume, flow_in, flow_out, dt] = ...
    wis_data("20210126_step_gate3_4_s255_no_intake.csv", wis);

figure(3);
title("Water heights");
plot(timing/1000, water_heights);
legend("s1", "s2", "s3", "s4", "s5", "s6", "s7")



%%

% data for identification

figure(10);
[input1a, input1b] = get_flows_for_pool(flow_in, flow_out, 3, true);
input1b = flow_out(:,3);

output1 = water_heights(:, 7) - water_heights(1, 7);


plot(input1a)
hold on;
plot(input1b)

[M, ~] = size(input1a);
input1 = zeros(M,1);
input1(1:111) = linspace(0,0.003228, 111);
input1(112:405) = linspace(0.003228, 0, 405 - 112 + 1);

output1a = output1;
% fake_output
output1(1) = 0.05; % start at
for i=2:M
    output1(i) = output1(i-1) + (input1(i-1) * dt) / wis.area3;
end


% figure(20)
% plot(output1)

%plot(input1);

%%
% plot(sqrt(water_heights(:,5) - water_heights(:,6)))
% hold on;
% plot(input1 * 200)
% 
% 
% input1 = height_to_flow(water_heights(:,5), water_heights(:,6));

%%

dt1 = dt;

[water_heights, timing, delta_volume, flow_in, flow_out, dt] = ...
    wis_data("20210126_step_gate3_4_s100_no_intake.csv", wis);


input2a = flow_in(:,4);
input2b = flow_out(:,3);

output2 = water_heights(:, 7) - water_heights(1, 7);
figure(11)

plot(input2a)
hold on;
plot(input2b)

[M, ~] = size(input2a);
input2 = zeros(M,1);
input2(1:89) = linspace(0,0.00186, 89);
input2(90:652) = linspace(0.00186, 0, 652 - 90 + 1);
% input2(1:19) = linspace(0,0.0006897, 19);
% input2(20:2136) = linspace(0.0006897, 0, 2136 - 20 + 1);

plot(input2);

output2a = output2;
% fake_output
output2(1) = 0.05; % start at
for i=2:M
    output2(i) = output2(i-1) + (input2(i-1) * dt) / wis.area3;
end
% 
% figure(20)
% plot(output1)

%%

dt2 = dt;


% TODO: check negative values + scale values
%input2 = height_to_flow(water_heights(:,5), water_heights(:,6));

% shift to compensate for delay

[input1a, output1a] = shift_delay(input1a(145:350), output1a(145:350), 20);
[input2a, output2a] = shift_delay(input2a(100:590), output2a(100:590), 20);

output1a = output1a - output1a(1);
output2a = output2a - output2a(1);


ze = iddata(output1a,input1a*1000,dt1); % experimemnt
figure(4)
plot(ze)

%%

zv = iddata(output2a,input2a*1000,dt2); % validate
figure(12)
plot(zv)

%%

Ge = spa(ze);
figure(5)
bode(Ge)

%%
% model estimation
Mimp = impulseest(ze,60); % 3rd order instead of 60?
figure(6)
% step response
step(Mimp)

%%


Mimp2 = impulseest(zv,60); % 3rd order instead of 60?
figure(13)
% step response
step(Mimp2)
%%

disp('delay')
delayest(ze)
delayest(zv)

%%

% Estimating delay not working correctly on this data
% Delays read from step graphs (response time at sensor on other side after 
% opening the gate) => rough estimates!

% pool1 1.2s (20 samples)
% pool2 0.4s (7 samples)
% pool3 1.4s (22 samples)

Opt = tfestOptions('Display','on');
np = [3]; % 3rd order
nz = [2]
ioDelay = [0];



mtf = tfest(zv,np,nz,ioDelay,Opt);

%%
figure(14);
compare(ze,mtf)

% 
% clf
% mi = impulseest(pool); % non-parametric (FIR) model
% showConfidence(impulseplot(mi),3); %impulse response with 3 standard
%                                    %deviations confidence region
%                                    
%                                    
%           