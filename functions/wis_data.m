function [data] = wis_data(csv_file, wis, pool, type, description, dt)
% load and augment csv data of the WIS lab setup

    % check if descriptive arguments have been set
    if nargin < 3
        pool = 0;
        type = "no type";
        description = "no description";
    end

    % load data
    pool_data = readmatrix(sprintf("data/%s", csv_file));
    
%     % resample data
%     r = 2;
%     y = zeros(ceil(size(pool_data,1)/r), size(pool_data,2));
%     for i = 1:size(pool_data,2)
%         x = pool_data(:, i); 
%         y(:, i) = decimate(x,r); 
%     end
%     pool_data = y;

    % convert pressure sensor data to water level 
    water_levels = pool_data(:, [3,4,7,8,11,12, 15]) .* wis.a + wis.b; % cm
    water_levels = water_levels / 100; % m
    
    
    
    % extract timing
    timing = pool_data(:, 1); % ms 
    
    
    [M, ~] = size(water_levels);

    % Manually set dt to avoid rounding errors
    if nargin < 6
        % use average sample time as timestep
        dt = ((timing(M) - timing(1)) / (M-1)) / 1000; % s
    end
    
    % init matrices to store volume change, flows and height difference
    delta_volume = zeros(M,4); % m^3
    flow_in = zeros(M,4); % m^3/sec
    flow_out = zeros(M,4); % m^3/sec
    delta_height = zeros(M,4); % m

    % TODO: can this be vectorized?
    % TODO: write more compact when done (for now this is easier to debug)
    for i = 2:M
        % calculate change in height per pool (use average for pools with 2
        % sensors)
        dh0 = water_levels(i, 1) - water_levels(i-1, 1); 
        dh1 = (water_levels(i, 2) - water_levels(i-1, 2) + water_levels(i, 3) - water_levels(i-1, 3)) / 2; 
        dh2 = (water_levels(i, 4) - water_levels(i-1, 4) + water_levels(i, 5) - water_levels(i-1, 5)) / 2; 
        dh3 = (water_levels(i, 6) - water_levels(i-1, 6) + water_levels(i, 7) - water_levels(i-1, 7)) / 2; 

        dv0 = dh0 * wis.area0;
        dv1 = dh1 * wis.area1;
        dv2 = dh2 * wis.area2;
        dv3 = dh3 * wis.area3;

        delta_volume(i,:) = [dv0 dv1 dv2 dv3];

        % calculate flow that would have lead to this change in height
        % (ZOH) based on pools before fore the flow_in and based on the
        % pools after this one for the flow_out
        
        fi0 = 0; % unknown
        fi1 = -dv0 / dt + fi0;
        fi2 = -dv1 / dt + fi1;
        fi3 = -dv2 / dt + fi2;

        flow_in(i-1,:) = [fi0 fi1 fi2 fi3];

        fo3 = 0; % unknown
        fo2 = dv3 / dt + fo3;
        fo1 = dv2 / dt + fo2;
        fo0 = dv1 / dt + fo1;

        flow_out(i-1,:) = [fo0 fo1 fo2 fo3];
    end
    
    % TODO: can this be vectorized?
    
    for i = 1:M
        % calculate height difference between pools
        delta_height(i,1) = water_levels(i, 1) - water_levels(i, 2);
        delta_height(i,2) = water_levels(i, 3) - water_levels(i, 4);
        delta_height(i,3) = water_levels(i, 5) - water_levels(i, 6);
        % TODO: this should be the 'head over the gate' from Cantoni
        % because gate 4 is an overshot gate
        delta_height(i,4) = water_levels(i, 7) - 0; % = gate height
    end
    
    % TODO: tune speed of the servos
    % apply rate limiter to actuator values
    actuators = pool_data(:, [5,9,13,17]);
    
    steps_per_sec = 255/8;
    steps_per_sample = steps_per_sec * dt;
    
    for i = 2:M
        % TODO: vectorize?
        for j = 1:4
            if actuators(i, j) >  actuators(i-1, j) + steps_per_sample
                actuators(i, j) = actuators(i-1, j) + steps_per_sample;
            end
            % TODO: also add closing
        end
    end
    
    data.water_levels = water_levels;
    data.timing = timing;
    data.delta_volume = delta_volume;
    data.delta_height = delta_height;
    data.flow_in = flow_in;
    data.flow_out = flow_out;
    data.dt = dt;
    data.sensors = pool_data(:, [3,4,7,8,11,12,15]);
    data.actuators = actuators;
    data.filename = csv_file;
    data.pool = pool;
    data.type = type;
    data.description = description;
end

