clear all;

calibration_data = readmatrix("20201105_fill_all_5_33cm.csv");

figure(1);
plot(calibration_data(:, 2:8));

[M, N] = size(calibration_data);

lowest_point = 5;
highest_point = 33;

cvx_begin quiet

variable a(7)  
variable b(7)  

objective = 0;

% 3 has high accuracy and is near source, so use as reference
a(2) * calibration_data(1, 3) + b(2) == lowest_point;
a(2) * calibration_data(2027, 3) + b(2) == highest_point;

% set constraints on the different type of sensors
a(1) == a(3);
a(1) == a(4);
a(1) == a(5);

a(2) == a(6);
a(2) == a(7);


% using all points is a bit overkill, so skip to speed up
for i = 1:100:M
    i
    for j = 1:6
        x1 = a(j) * calibration_data(i, j+1) + b(j);
        x2 = a(j+1) * calibration_data(i, j+2) + b(j+1);
        error = (x1-x2)^2;
        objective = objective + error;
    end
end

minimize(objective);

cvx_end

%% plot again after calibration

figure(2);
calibrated_sensors = calibration_data(:, 2:8) .* a' + b';

plot(calibrated_sensors);


%% convert identification data to mm

identification_data = readmatrix("20201105_equilibrium_gates_1cm.csv");

figure(3);
calibrated_identification = identification_data(:, 2:8) .* a' + b';

plot(calibrated_identification)

%% add flows through gate 1, 2 and 3
area1 = 0.1853; %m2
area2 = 0.1187; %m2
area3 = 0.2279; %m2

dt = 0.5; %sec

% TODO: can this be vectorized?
[M, N] = size(calibrated_identification);
calibrated_identification = [calibrated_identification zeros(M, 3)];

for i = 2:M
    % use average in pool to estimate height difference between timesteps
    dh3 = (((calibrated_identification(i, 6) + calibrated_identification(i, 7)) / 2) - ...
        ((calibrated_identification(i-1, 6) + calibrated_identification(i-1, 7)) / 2)) * 0.001; %m
    dh2 = (((calibrated_identification(i, 4) + calibrated_identification(i, 5)) / 2) - ...
        ((calibrated_identification(i-1, 4) + calibrated_identification(i-1, 5)) / 2)) * 0.001; %m
    dh1 = (((calibrated_identification(i, 2) + calibrated_identification(i, 3)) / 2) - ...
        ((calibrated_identification(i-1, 2) + calibrated_identification(i-1, 3)) / 2)) * 0.001; %m
    
    % calculate flow over gate from changed water level in next gates
    gate3 = (dh3 * area3) / dt; %m3/sec
    gate2 = (dh2 * area2) / dt + gate3; %m3/sec
    gate1 = (dh1 * area1) / dt + gate2 + gate3; %m3/sec
    
    % update data
   calibrated_identification(i, 8:10) = [gate1 gate2 gate3]; 
end

plot(lowpass(calibrated_identification(:, 8:10), 0.01, 2))


