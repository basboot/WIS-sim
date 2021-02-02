function [data] = wis_data(csv_file, wis)
% load and augment csv data of the WIS lab setup

    pool_data = readmatrix(csv_file);

    water_heights = pool_data(:, [3,4,7,8,11,12, 15]) .* wis.a + wis.b; % cm
    water_heights = water_heights / 100; % m
    timing = pool_data(:, 1);

%     water_heights_s = water_heights;
%     for I = 1:7
%         water_heights_s(:,I) = smooth(water_heights(:,I), 100);
%     end

    % TODO: no smoothing anymore, cleanup
    water_heights_s = water_heights;
    
    
    % calculate flows
    [M, ~] = size(water_heights);

    % average sample time
    dt = ((timing(M) - timing(1)) / (M-1)) / 1000;

    %  volume change and flows
    delta_volume = zeros(M,4); % m3
    flow_in = zeros(M,4); % m3/sec
    flow_out = zeros(M,4); % m3/sec
    delta_height = zeros(M,4);

    % TODO: can this be vectorized?
    % TODO: write more compact when done (for now this is easier to debug)
    for i = 2:M
        % calculate change in height per pool (use average for pools with 2
        % sensors)
        dh0 = water_heights_s(i, 1) - water_heights_s(i-1, 1); 
        dh1 = (water_heights_s(i, 2) - water_heights_s(i-1, 2) + water_heights_s(i, 3) - water_heights_s(i-1, 3)) / 2; 
        dh2 = (water_heights_s(i, 4) - water_heights_s(i-1, 4) + water_heights_s(i, 5) - water_heights_s(i-1, 5)) / 2; 
        dh3 = (water_heights_s(i, 6) - water_heights_s(i-1, 6) + water_heights_s(i, 7) - water_heights_s(i-1, 7)) / 2; 

        dv0 = dh0 * wis.area0;
        dv1 = dh1 * wis.area1;
        dv2 = dh2 * wis.area2;
        dv3 = dh3 * wis.area3;

        delta_volume(i,:) = [dv0 dv1 dv2 dv3];

        fi0 = 0; % unknown
        fi1 = -dv0 / dt + fi0;
        fi2 = -dv1 / dt + fi1;
        fi3 = -dv2 / dt + fi2;

        flow_in(i,:) = [fi0 fi1 fi2 fi3];

        fo3 = 0; % unknown
        fo2 = dv3 / dt + fo3;
        fo1 = dv2 / dt + fo2;
        fo0 = dv1 / dt + fo1;

        flow_out(i,:) = [fo0 fo1 fo2 fo3];
        
        %water_heights = water_heights_s;
    end
    
    % TODO: can this be vectorized?
    
    for i = 1:M
        % calculate height difference between pools
        delta_height(i,1) = water_heights_s(i, 1) - water_heights_s(i, 2);
        delta_height(i,2) = water_heights_s(i, 3) - water_heights_s(i, 4);
        delta_height(i,3) = water_heights_s(i, 5) - water_heights_s(i, 6);
        % TODO: this should be the overshot from Cantoni
        delta_height(i,4) = water_heights_s(i, 7) - 0; % = gate height
    end
    
    data.water_heights = water_heights;
    data.timing = timing;
    data.delta_volume = delta_volume;
    data.delta_height = delta_height;
    data.flow_in = flow_in;
    data.flow_out = flow_out;
    data.dt = dt;
    data.sensors = pool_data(:, [3,4,7,8,11,12,15]);
    data.actuators = pool_data(:, [5,9,13,17]);

end

