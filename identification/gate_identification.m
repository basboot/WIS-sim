%% identify gate parameters

% simplified model:
% flow = sign(dh) * gate_opening * K * sqrt(abs(dh))
% K not constant, probably due to leakage so a better approximation will
% be:
% flow = sign(dh) * gate_opening * K(gate_opening) * sqrt(abs(dh))


% Identify 3 gates
% Identify 3 settings
%   - use only part where gate is fully open
% Plot results

% TODO: automate data selection

gate(1).data = [pool_data(1), pool_data(2), pool_data(3)];
gate(1).setting = [25, 100, 255];
gate(1).flow_to_use = [2, 3, 3];
gate(1).k = zeros(1, size(gate(1).setting, 2));

gate(2).data = [pool_data(4), pool_data(5), pool_data(6)];
gate(2).setting = [25, 100, 255];
gate(2).flow_to_use = [2, 3, 3];
gate(2).k = zeros(1, size(gate(2).setting, 2));

gate(3).data = [pool_data(7), pool_data(8), pool_data(9)];
gate(3).setting = [25, 100, 255];
gate(3).flow_to_use = [2, 3, 3];
gate(3).k = zeros(1, size(gate(3).setting, 2));

for gate_number = 1: 3
    for i = 1:size(gate(gate_number).setting, 2);


        gate(gate_number).k(i) = calculate_flow_constant(gate(gate_number).data(i), gate_number, gate(gate_number).setting(i), true, gate(gate_number).flow_to_use(i), false);
    end
    
end


%% Plot results
figure();
hold on;
for gate_number = 1: 3
    scatter(gate(gate_number).setting, gate(gate_number).k);   
end


coefficients = polyfit([gate(1).setting gate(2).setting gate(3).setting],...
    [gate(1).k gate(2).k gate(3).k], 1);

a_gate = coefficients (1);
b_gate = coefficients (2);

plot([0 255], [0 255] .* a_gate + b_gate);

legend("gate1", "gate2", "gate3", "linear approximation");
xlabel("servo setting");
ylabel("flow constant");

% store results in wis data
wis.a_gate = a_gate;
wis.b_gate = b_gate;



