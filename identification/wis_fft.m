%% NOT FINSHED, AND NOT USED, YET

%% Frequency analysis of the pools


%% Load experiment data
gate3 = ...
    wis_data("20210202_step_gate3_4_s255_no_intake2.csv", wis);

% show time plot
figure(1);
plot(gate3.timing/1000, gate3.water_levels);
title("Water levels");
legend("s1", "s2", "s3", "s4", "s5", "s6", "s7")
xlabel("time [s]");
ylabel("water level [m]");


Fs = 1/gate3.dt;        % Sampling frequency                    
T = gate3.dt;       % Sampling period       
L = size(gate3.water_levels,1);    % Length of signal
t = (0:L-1)*T;        % Time vector
X = gate3.water_levels(:, 7)';

figure(1);
plot(1000*t(1:50),X(1:50))
title('Signal')
xlabel('t (milliseconds)')
ylabel('X(t)')

Y = fft(X);

P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

figure(2);
f = Fs*(0:(L/2))/L;
plot(f,P1) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')