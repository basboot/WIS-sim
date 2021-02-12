function [k] = calculate_flow_constant(gate_data, gate_to_identify, gate_setting, skip_opening, flow_to_use, showPlot)
% flow to use
% 0 = none, 1 = in, 2 = out, 3 = both
assert(flow_to_use > 0, "WARNING: in or outflow must be used to identify the gate.");

% simplified model:
% flow = sign(dh) * gate_opening * K * sqrt(abs(dh))
% K not constant, probably due to leakage so a better approximation will
% be:
% flow = sign(dh) * gate_opening * K(gate_opening) * sqrt(abs(dh))

% show message as feedback that the script is still doing something
fprintf("Calculating flow constant for gate %d with setting %d.\n", gate_to_identify, gate_setting);

first_measurement_to_use = 1;
% Remove opening of the gate
if skip_opening
    first_measurement_to_use = find(gate_data.actuators(:, gate_to_identify) == gate_setting, 1, 'first');
end

last_measurement_to_use = size(gate_data.actuators(:, gate_to_identify),1)-500;

%%

% show time plot
if showPlot
    figure();
    plot(gate_data.timing/1000, gate_data.water_levels);
    title("Water levels");
    legend("s1", "s2", "s3", "s4", "s5", "s6", "s7")
    xlabel("time [s]");
    ylabel("water level [m]");
end

[f_in, f_out] = get_flows_for_pool(gate_data.flow_in, gate_data.flow_out, gate_to_identify, false);

%% show flows calculated from derivative
if showPlot
    figure(2);
    plot(f_in);
    hold on;
    plot(f_out);
    title("Flows calculated using the water levels")
    legend("q_{in}", "q_{out}");
end


%% show filtered flows, and flow estimated from difference in water height
% Filtered flows are for illustrative purposes only

if showPlot
    figure();
    plot(lowpass(f_in,1,1/gate_data.dt))
    hold on
    plot(lowpass(f_out,1,1/gate_data.dt))
end

% % TODO: use optimization tool for accurate estimation of K
% K = 0.00004;
% % 255 => 200
% 
% plot(gate_setting * K * sign(gate_data.delta_height(:,gate_to_identify)) .* sqrt(abs(gate_data.delta_height(:,gate_to_identify))))

cvx_clear

cvx_begin quiet

variable k 

objective = 0;

% optimize sum squared error between model and estimated flow
% use inflow for identification (1 or 3)
if bitand(flow_to_use, 1)
    objective = objective + sum((k * gate_data.actuators(first_measurement_to_use:last_measurement_to_use,gate_to_identify) .* sign(gate_data.delta_height(first_measurement_to_use:last_measurement_to_use,gate_to_identify)) .* sqrt(abs(gate_data.delta_height(first_measurement_to_use:last_measurement_to_use,gate_to_identify))) - f_in(first_measurement_to_use:last_measurement_to_use)).^2);
end

% use outflow for identification (2 or 3)
if bitand(flow_to_use, 2)
    objective = objective + sum((k * gate_data.actuators(first_measurement_to_use:last_measurement_to_use,gate_to_identify) .* sign(gate_data.delta_height(first_measurement_to_use:last_measurement_to_use,gate_to_identify)) .* sqrt(abs(gate_data.delta_height(first_measurement_to_use:last_measurement_to_use,gate_to_identify))) - f_out(first_measurement_to_use:last_measurement_to_use)).^2);
end

minimize(objective);

cvx_end

if showPlot
    plot(k * gate_data.actuators(:,gate_to_identify) .* sign(gate_data.delta_height(:,gate_to_identify)) .* sqrt(abs(gate_data.delta_height(:,gate_to_identify))));

    title("Filtered flows and model comparision")
    legend("q_{in}", "q_{out}", "q_{est}");
end

end

