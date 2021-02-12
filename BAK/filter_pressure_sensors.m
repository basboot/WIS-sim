dt = 0.1; % measurement speed

filter_data = readmatrix("20201201_sensors202_fill_gate100_bw0_5.csv");

plot(filter_data(:, 3:6))


%%

%% TODO: filter data before processing

[m_i, n_i] = size(filter_data);

filtered_data = zeros(m_i, n_i);


% butter
fc = 0.25; % cut off 
fs = 10; % sampling rate 

[b_b,b_a] = butter(1,fc/(fs/2)) %1st order

%freqz(b_b,b_a)


% process
filtered_data(:, 3:4) = filter(b_b,b_a,filter_data(:, 3:4));

%filtered_data(:, 3:4) = lowpass(filter_data(:, 3:4), 1, 10);

% %start averaging
% last_value = filter_data(1, 3:4);
% 
% w_last = 0.75;
% w_current = 1 - w_last;
% 
% for i = 1:m_i
%     current_value = filter_data(i, 3:4);
%     new_value = w_current * current_value + w_last * last_value;
%     filtered_data(i, 3:4) = new_value;
%     last_value = new_value;
% end

%start bw
x0 = zeros(1,2);
y0 = zeros(1,2);
x1 = zeros(1,2);
y1 = zeros(1,2);

% TODO: calculate correct values
%GAIN = 245.8523972;
%GAIN = 4.077683537;
GAIN = 0.0703; %b_b(2); %4.077683537; %1/0.2452; %4; % 1/B1/2

%A = 0.9918650376 
%A = 0.5095254495;
A = 0.8594; %-b_a(2); %0.5095254495; %0.5; % A2

% last_value = filter_data(1, 3:4); 
% 
% for i = 1:m_i
%     
%     current_value = filter_data(i, 3:4);
%     
%     x0 = x1; 
%     x1 = current_value * GAIN;
%     y0 = y1; 
%     y1 =   (x0 + x1) + (  A * y0);
%        
%     
%     filtered_data(i, 3:4) = y1;
% end


%%

figure(3);
calibrated_filter_data = filter_data(:, 3:4) .* a(3:4)' + b(3:4)';

calibrated_filtered_data = filtered_data(:, 3:4) .* a(3:4)' + b(3:4)';


plot(calibrated_filter_data)
hold on;
plot(calibrated_filtered_data)
title("filter(ed) data - time/water heigts");

