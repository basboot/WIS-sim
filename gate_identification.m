%% identify gate parameters

% simplified model:
% flow = sign(dh) * gate_opening * K * sqrt(abs(dh))

%% Load experiment data
gate3 = ...
    wis_data("20210202_step_gate3_4_s25_no_intake2.csv", wis);

% show time plot
figure(1);
plot(gate3.timing/1000, gate3.water_levels);
title("Water levels");
legend("s1", "s2", "s3", "s4", "s5", "s6", "s7")
xlabel("time [s]");
ylabel("water level [m]");

% % Not needed anymore => automated in wis_data
% %% change gate
% secs_to_open = 8;
% sps = 128;
% value_to_open = 255;
gate_to_open = 3;
% 
% gate3.actuators(1:secs_to_open * sps, gate_to_open) = linspace(0, value_to_open, secs_to_open * sps);

%% show flows calculated from derivative
figure(2);
[f_in, f_out] = get_flows_for_pool(gate3.flow_in, gate3.flow_out, 3, false);


plot(f_in);
hold on;
plot(f_out);
title("Flows calculated using the water levels")
legend("q_{in}", "q_{out}");

%% show filtered flows, and flow estimated from difference in water height
figure(3);
plot(lowpass(f_in,1,1/gate3.dt))
hold on
plot(lowpass(f_out,1,1/gate3.dt))

% TODO: use optimization tool for accurate estimation of K
K = 0.00004;
% 255 => 200
gate_opening = 200; 
plot(gate_opening * K * sign(gate3.delta_height(:,3)) .* sqrt(abs(gate3.delta_height(:,3))))

cvx_begin quiet

variable k 

objective = 0;

% sum squared error model and estimated inflow
objective = sum((k * gate3.actuators(:,gate_to_open) .* sign(gate3.delta_height(:,3)) .* sqrt(abs(gate3.delta_height(:,3))) - f_in).^2);
% sum squared error model and estimated outflow
%objective = objective + sum((k * gate3.actuators(:,gate_to_open) .* sign(gate3.delta_height(:,3)) .* sqrt(abs(gate3.delta_height(:,3))) - f_out).^2);

minimize(objective);

cvx_end

plot(k * gate3.actuators(:,gate_to_open) .* sign(gate3.delta_height(:,3)) .* sqrt(abs(gate3.delta_height(:,3))))


title("Filtered flows and model comparision")
legend("q_{in}", "q_{out}", "q_{est}", "q_{solv}");

%% TODO: automate...

% TODO: can this be vectorized?
coefficients = polyfit([25, 255], [5.0, 3.5], 1);
a_gate = coefficients (1);
b_gate = coefficients (2);
