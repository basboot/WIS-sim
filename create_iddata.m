function [ze] = create_iddata(data_experiment, showPlot)
% set dt manually to avoid rounding errors
%   Detailed explanation goes here

pool = data_experiment.pool;


if showPlot
    figure(1);
    title("Water levels");
    plot(data_experiment.timing/1000, data_experiment.water_levels);
    legend("s1", "s2", "s3", "s4", "s5", "s6", "s7")
    xlabel("time [s]");
    ylabel("water level [m]");
end

if showPlot
    figure(2);
end
[input1, ~] = get_flows_for_pool(data_experiment.flow_in, data_experiment.flow_out, pool, showPlot);

%f_est_exp = 200 * 0.00004 * sign(data_experiment.delta_height(:,3)) .* sqrt(abs(data_experiment.delta_height(:,3)));

% calculate using gate ident data
% use k to compensate for gate opening irregularities
% gate_to_open = 3;
% k = (data_experiment.actuators(:,gate_to_open) .* a_gate + b_gate) * 10^-5;
% 
% f_est_exp = k .* data_experiment.actuators(:,gate_to_open) .* sign(data_experiment.delta_height(:,3)) .* sqrt(abs(data_experiment.delta_height(:,3)));

% if showPlot
%     plot(f_est_exp)
% end
%input1 = f_est_exp;

output1 = data_experiment.water_levels(:, pool*2+1);
dt1 = data_experiment.dt;

%% shift data (for delay) and select most reliable part

% Estimating delay not working correctly on this data
% Delays read from step graphs (response time at sensor on other side after 
% opening the gate) => rough estimates!

% pool1 1.2s (20 samples) 16sps
% pool2 0.4s (7 samples)
% pool3 1.4s (22 samples)
% 
% [input1, output1] = shift_delay(input1(119:386), output1(119:386), 0);
% [input2, output2] = shift_delay(input2(119:400), output2(119:400), 0);



% crop
input1 = input1(1:end);
output1 = output1(1:end);


% remove 'datum'
output1 = output1 - output1(1);

%% experiment

% input converted to l/sec for better scaling!
ze = iddata(output1,input1*1000,dt1); 
if showPlot
    figure(4)
    plot(ze)
end

% %% bode experiment
% 
% Ge = spa(ze);
% figure(5)
% bode(Ge)

%%
% model estimation, experiment
Mimp = impulseest(ze,60); 
if showPlot
    figure(6)
    % step response
    step(Mimp)
end

%% estimate delays
disp('delay')
delayest(ze)

ze.ExperimentName = sprintf("%s-%s", data_experiment.description, data_experiment.type);
ze.Name = sprintf("%s-%s", data_experiment.filename, data_experiment.type);


end

