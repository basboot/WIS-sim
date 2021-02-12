%% calibration.m

% Calibrate the presure sensors (needed to convert pressure to water level)

CalibrationData = readmatrix("../data/20210202_calib_5cm.csv");

% figure(1);
% plot(calibration_data(:, [3,4,7,8,11,12, 15]));
% title("calibration data raw (5cm)");

low_value = mean(CalibrationData(:, [3,4,7,8,11,12, 15]));
low = 5;

CalibrationData = readmatrix("../data/20210202_calib_31cm.csv");

% figure(2);
% plot(calibration_data(:, [3,4,7,8,11,12, 15]));
% title("calibration data raw (31cm)");

high_value = mean(CalibrationData(:, [3,4,7,8,11,12, 15]));
high = 31;

% fit data

% TODO: can this be vectorized?
[M, N] = size(low_value);
a = zeros(1,N);
b = zeros(1,N);
for i = 1:N
    coefficients = polyfit([low_value(i), high_value(i)], [low, high], 1);
    a(i) = coefficients (1);
    b(i) = coefficients (2);
end

Wis.a = a;
Wis.b = b;

%%

% [M, N] = size(calibration_data);
% 
% % first
% lowest_point = 33;
% % last
% highest_point = 5; 
% 
% cvx_begin quiet
% 
% variable a(7)  
% variable b(7)  
% 
% objective = 0;
% 
% % 3 has high accuracy and is near source, so use as reference
% a(2) * calibration_data(1, 3) + b(2) == lowest_point;
% a(2) * calibration_data(M-1, 3) + b(2) == highest_point;
% 
% % set constraints on the different type of sensors
% a(1) == a(3);
% a(1) == a(4);
% a(1) == a(5);
% 
% a(2) == a(6);
% a(2) == a(7);
% 
% 
% % using all points is a bit overkill, so skip to speed up
% for i = 1:100:M
%     i
%     for j = 1:6
%         x1 = a(j) * calibration_data(i, j+1) + b(j);
%         x2 = a(j+1) * calibration_data(i, j+2) + b(j+1);
%         error = (x1-x2)^2;
%         objective = objective + error;
%     end
% end
% 
% minimize(objective);
% 
% cvx_end
% 
% %% plot again after calibration
% 
% figure(2);
% calibrated_sensors = calibration_data(:, 2:8) .* a' + b';
% 
% plot(calibrated_sensors);
% title("calibration data - calibrated");