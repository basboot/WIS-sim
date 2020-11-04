clear all;

calibration_data = readmatrix("20201103_fill_all.csv");

figure(1);
plot(calibration_data(:, 2:8))

[M, N] = size(calibration_data)

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

%% plot again afyet calibration

figure(2);
sensors = calibration_data(:, 2:8) .* a' + b'

plot(sensors)

%% 


% a =
% 
%     0.0183
%     0.0053
%     0.0183
%     0.0183
%     0.0183
%     0.0053
%     0.0053
%     
% b =
% 
%  -183.3012
%   -64.7004
%  -183.7719
%  -181.5880
%  -179.6035
%   -58.2928
%   -58.4973
