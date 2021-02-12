%% gate_identification.m

% Identify gate parameters

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

Gate(1).data = [PoolData(1), PoolData(2), PoolData(3)];
Gate(1).setting = [25, 100, 255];
Gate(1).flow_to_use = [2, 3, 3];
Gate(1).k = zeros(1, size(Gate(1).setting, 2));

Gate(2).data = [PoolData(4), PoolData(5), PoolData(6)];
Gate(2).setting = [25, 100, 255];
Gate(2).flow_to_use = [2, 3, 3];
Gate(2).k = zeros(1, size(Gate(2).setting, 2));

Gate(3).data = [PoolData(7), PoolData(8), PoolData(9)];
Gate(3).setting = [25, 100, 255];
Gate(3).flow_to_use = [2, 3, 3];
Gate(3).k = zeros(1, size(Gate(3).setting, 2));

for gate_number = 1: 3
    for i = 1:size(Gate(gate_number).setting, 2)
        Gate(gate_number).k(i) = calculateFlowConstant(Gate(gate_number).data(i), gate_number, Gate(gate_number).setting(i), true, Gate(gate_number).flow_to_use(i), false);
    end
    
end


%% Plot results
figure();
hold on;
for gate_number = 1: 3
    scatter(Gate(gate_number).setting, Gate(gate_number).k);   
end


coefficients = polyfit([Gate(1).setting Gate(2).setting Gate(3).setting],...
    [Gate(1).k Gate(2).k Gate(3).k], 1);

a_gate = coefficients (1);
b_gate = coefficients (2);

plot([0 255], [0 255] .* a_gate + b_gate);

legend("gate1", "gate2", "gate3", "linear approximation");
xlabel("servo setting");
ylabel("flow constant");

% store results in wis data
Wis.a_gate = a_gate;
Wis.b_gate = b_gate;



