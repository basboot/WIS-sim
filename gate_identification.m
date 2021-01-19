dt = 0.5; % sec
servo_speed = 32; % steps/sec

%% convert identification data to mm

%identification_data = readmatrix("20201117_gate2faster_10-250.csv");
identification_data = readmatrix("20201201_fill_all.csv");

%% TODO: filter data before processing

% butter
fc = 0.1; % cut off 
fs = 10; % sampling rate 

[b_b,b_a] = butter(1,fc/(fs/2)); %1st order
%freqz(b,a)

% process
%identification_data(:, 2:8) = filter(b_b,b_a,identification_data(:, 2:8));

%identification_data(:, 2:8) = lowpass(identification_data(:, 2:8), 0.02, 2);

[m_i, n_i] = size(identification_data);

% %start
% last_value = identification_data(1, 2:8);
% 
% w_last = 0.99;
% w_current = 1 - w_last;
% 
% for i = 1:m_i
%     current_value = identification_data(i, 2:8);
%     new_value = w_current * current_value + w_last * last_value;
%     identification_data(i, 2:8) = new_value;
%     last_value = new_value;
% end


%%

figure(3);
calibrated_identification = identification_data(:, 2:8) .* a' + b';


plot(calibrated_identification)
title("identification data - time/water heigts");

%% add flows through gate 1, 2 and 3
area1 = 0.1853; %m2
area2 = 0.1187; %m2
area3 = 0.2279; %m2

[M, N] = size(calibrated_identification);

% add gate openings and make room for dh and flow over gate
calibrated_identification = [calibrated_identification zeros(M, 6) identification_data(:, 9:12)];

gates = identification_data(1, 9:12);

% TODO: can this be vectorized?
for i = 2:M
    % use average in pool to estimate height difference between two timesteps
    dh3 = (((calibrated_identification(i, 6) + calibrated_identification(i, 7)) / 2) - ...
        ((calibrated_identification(i-1, 6) + calibrated_identification(i-1, 7)) / 2)) * 0.001; %m
    dh2 = (((calibrated_identification(i, 4) + calibrated_identification(i, 5)) / 2) - ...
        ((calibrated_identification(i-1, 4) + calibrated_identification(i-1, 5)) / 2)) * 0.001; %m
    dh1 = (((calibrated_identification(i, 2) + calibrated_identification(i, 3)) / 2) - ...
        ((calibrated_identification(i-1, 2) + calibrated_identification(i-1, 3)) / 2)) * 0.001; %m
    
    % calculate flow over gate from changed water level in next gates
    gate3 = (dh3 * area3) / dt; %m3/sec
    gate2 = (dh2 * area2) / dt + gate3; %m3/sec
    gate1 = (dh1 * area1) / dt + gate2 + gate3; %m3/sec
    
    % calculate height difference on both sides of the gate (use only
    % sensor near the gate
    dg1 = calibrated_identification(i, 1) - calibrated_identification(i, 2);
    dg2 = calibrated_identification(i, 3) - calibrated_identification(i, 4);
    dg3 = calibrated_identification(i, 5) - calibrated_identification(i, 6);
    
    % update data
    calibrated_identification(i, 8:13) = [dg1 dg2 dg3 gate1 gate2 gate3]; 
   
    %% TODO
    % Process data to fix gate opening speed (do not fix closing
    % yet!)
    % check is opening is larger than speed allowed
    % TODO: vectorize?
    for j = 1:4
        if calibrated_identification(i, 13+j) > gates(j) + dt * servo_speed
            calibrated_identification(i, 13+j) = gates(j) + dt * servo_speed;
        end
    end
    
    % update current
    gates = calibrated_identification(i, 14:17);


end

%%
% Delete rows with zero gate for a fixed column
ids = calibrated_identification(:,15) == 0;
calibrated_identification(ids, :) = [];

%%
figure(4);
plot((calibrated_identification(:, 9)))
title("identification data - time/heigt difference");

figure(5);
plot(calibrated_identification(:, 12))
title("identification data - time/flow");

%% flow vs height
[m_o, n_o] = size(calibrated_identification) % TODO: think of a better solution to cleanup data
overflow = m_o; % sample at which overflow at gate 4 starts

figure(6);
% flow_vs_height1 = [calibrated_identification(1:overflow, 8), calibrated_identification(1:overflow, 11)];
% % sort by height for plotting
% %[~,idx] = sort(flow_vs_height1(:,1)); 
% %flow_vs_height1 = flow_vs_height1(idx,:);
% 
% scatter(flow_vs_height1(:,1), flow_vs_height1(:,2), 'filled')

flow_vs_height2 = [calibrated_identification(1:overflow, 9), calibrated_identification(1:overflow, 12)];
% sort by height for plotting
%[~,idx] = sort(flow_vs_height2(:,1)); 
%flow_vs_height2 = flow_vs_height2(idx,:);
title("identification data - height difference/flow");

hold on;

scatter(flow_vs_height2(:,1), flow_vs_height2(:,2), 'filled')

% hold on;
% flow_vs_height3 = [calibrated_identification(1:overflow, 10), calibrated_identification(1:overflow, 13)];
% % sort by height for plotting
% %[~,idx] = sort(flow_vs_height3(:,1)); 
% %flow_vs_height3 = flow_vs_height3(idx,:);
% 
% scatter(flow_vs_height3(:,1), flow_vs_height3(:,2), 'filled')

%%
figure(7)
flow_vs_height_opening2 = flow_vs_height2; % just a copy for the first column

% TODO: use lowpass filter? or remove this (and use all data!!!)
%flow_vs_height_opening2(:,2) = flow_vs_height2(:,2) ./ (calibrated_identification(1:overflow, 15) + 15);

% filter before splice
flow_vs_height_opening2(:,2) = flow_vs_height_opening2(:,2);

% splice
flow_vs_height_opening2_10 = flow_vs_height_opening2(20:381,:);

scatter(flow_vs_height_opening2_10(:,1), flow_vs_height_opening2_10(:,2) ./ (10+10), 'filled')
hold on;

flow_vs_height_opening2_25 = flow_vs_height_opening2(395:565,:);

scatter(flow_vs_height_opening2_25(:,1), flow_vs_height_opening2_25(:,2) ./ (25+15), 'filled')

flow_vs_height_opening2_50 = flow_vs_height_opening2(580:708,:);

scatter(flow_vs_height_opening2_50(:,1), flow_vs_height_opening2_50(:,2) ./ (50+15), 'filled')

flow_vs_height_opening2_100 = flow_vs_height_opening2(724:773,:);

scatter(flow_vs_height_opening2_100(:,1), flow_vs_height_opening2_100(:,2) ./ (100+15), 'filled')

flow_vs_height_opening2_150 = flow_vs_height_opening2(790:821,:);

scatter(flow_vs_height_opening2_150(:,1), flow_vs_height_opening2_150(:,2) ./ (150+15), 'filled')

flow_vs_height_opening2_200 = flow_vs_height_opening2(839:858,:);

scatter(flow_vs_height_opening2_200(:,1), flow_vs_height_opening2_200(:,2) ./ (200+15), 'filled')

flow_vs_height_opening2_250 = flow_vs_height_opening2(876:895,:);

scatter(flow_vs_height_opening2_250(:,1), flow_vs_height_opening2_250(:,2) ./ (250+15) , 'filled')

% compare to square root by hand
heights = 0:30
scatter(heights, sqrt(heights) / 2750000, 'filled')

% max

%sqrt(heights) / 2750000 * 265 : assume max 60cm, and gate 300
disp('max flow');
sqrt(60) / 275000 * 300


% so formula will be:
% flow = corrected_gate_steps / 2750000 * sqrt(height) (m3/sec)


% Assume model in form K * gate_opening * sqrt(height)
% Ce unknown

 

% % TODO: problem not convex => use Yalmip
% 
% yalmip('clear')
% 
% K = sdpvar(1);
% 
% constraints = []
% 
% % gate 10
% [m_f, n_f] = size(flow_vs_height_opening2_10);
% objective = 0;
% for i = 1:m_f
%     if flow_vs_height_opening2_10(i,2) < 0
%         continue;
%     end
%     s_height = sqrt(flow_vs_height_opening2_10(i,2));
%     gate = 20;
%     y = K * gate * s_height;
%     error = (y-flow_vs_height_opening2_10(i,1))^2;
%     objective = objective + error;
% end
% 
% [m_f, n_f] = size(flow_vs_height_opening2_25);
% for i = 1:m_f
%     if flow_vs_height_opening2_25(i,2) < 0
%         continue;
%     end
%     s_height = sqrt(flow_vs_height_opening2_25(i,2));
%     gate = 40;
%     y = K * gate * s_height;
%     error = (y-flow_vs_height_opening2_25(i,1))^2;
%     objective = objective + error;
% end
% 
% optimize(constraints, objective, sdpsettings('debug',1, 'solver','mosek')); 
% 
% value(K)
% 
% 
