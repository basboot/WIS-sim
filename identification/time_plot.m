%% time_plot.m

% Show a graph to visualy check the calibration

PlotData = ...
    createWisData("20210302_step_gate1_2_s255_slow_no_intake.csv", ...
    Wis);

PlotData = ...
    createWisData("20210316_noise_test_pump_silent_servos.csv", ...
    Wis);

PlotData = ...
    createWisData("20210323_noise_test_pump_silent.csv", ...
    Wis);

% PlotData.timing = PlotData.timing(3000:end, :);
% PlotData.water_levels = PlotData.water_levels(3000:end, :);

figure();
title("Water levels");
plot(PlotData.timing/1000, PlotData.water_levels);
legend("s1", "s2", "s3", "s4", "s5", "s6", "s7")

% frequency analysis
Fs = 128;            % Sampling frequency                    
T = 1/Fs;             % Sampling period       
L = size(PlotData.water_levels, 1);             % Length of signal
t = (0:L-1)*T;        % Time vector

Y = fft(PlotData.water_levels);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = Fs*(0:(L/2))/L;
figure();
plot(f,P1) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')

% butter
fc = 1; % cut off 
fs = 128; % sampling rate 

[b_b,b_a] = butter(1,fc/(fs/2)); %1st order

%freqz(b_b,b_a)


% process data

% matlab filter
%filtered_levels = filter(b_b,b_a,PlotData.water_levels);

% c-style filter
[m_i, n_i] = size(PlotData.water_levels);
filtered_levels = zeros(m_i, n_i);

x0 = zeros(1,n_i);
y0 = zeros(1,n_i);
x1 = zeros(1,n_i);
y1 = zeros(1,n_i);

GAIN = 1/32-1/128; %0.0240; % 0.0240 b_b(2); %b_b(2)

A = 1-1/32-1/64; %0.9521; % 0.9521 -b_a(2); %-b_a(2)
 
for i = 1:m_i
    
    current_value = PlotData.water_levels(i, :);
    
    x0 = x1; 
    x1 = current_value * GAIN;
    y0 = y1; 
    y1 =   (x0 + x1) + (  A * y0);
       
    
    filtered_levels(i, :) = y1;
end


figure();
title("Water levels (filtered)");
plot(PlotData.timing(200:end,:)/1000, filtered_levels(200:end, :));
legend("s1", "s2", "s3", "s4", "s5", "s6", "s7")

% frequency analysis
Fs = 128;            % Sampling frequency                    
T = 1/Fs;             % Sampling period       
L = size(filtered_levels(200:end, :),1);             % Length of signal
t = (0:L-1)*T;        % Time vector

Y = fft(filtered_levels(200:end, :));
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = Fs*(0:(L/2))/L;
figure();
plot(f,P1) 
title('Single-Sided Amplitude Spectrum of X(t) after filtering')
xlabel('f (Hz)')
ylabel('|P1(f)|')

                  