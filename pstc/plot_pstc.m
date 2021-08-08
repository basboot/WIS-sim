%% 

% Time series from Simulink, with Data for 3 pools / FFs
% flow, level, servo, epoch, radio_on

target = [0.25 0.20 0.15];
showPlots = true;

useCachedData = true;

nExperiments = 1;
experiments{1} = ["20210806pstc_simulation", "PSTC", "(test)", "20210806pstc_controller"];
% filename sim, sigma, epsilon, filename contr

for iExperiments = 1: nExperiments

    experimentNameSim = experiments{iExperiments}(1);
    experimentNameCtrl = experiments{iExperiments}(4);
    
    fileName1 = sprintf('%s.mat', experimentNameSim);
    fileName2 = sprintf('%s.mat', experimentNameCtrl);

    
    %% Load cached data
    try
        load(fileName1);
    catch
        assert(false, "Data file sim does not exist");
    end


    if showPlots   
        figure();
        plot(y_log');
        hold on;

        xlabel('time (s)')
        ylabel('level (m)')
        legend('pool1', 'pool2', 'pool3');

        yline(0.25,'-','reference 1', 'LabelHorizontalAlignment', 'left', 'HandleVisibility','off');
        yline(0.20,'-','reference 2', 'LabelHorizontalAlignment', 'left', 'HandleVisibility','off');
        yline(0.15,'-','reference 3', 'LabelHorizontalAlignment', 'left', 'HandleVisibility','off');

        saveFigureEps(sprintf("%s-levels", experimentNameSim));

        title("Water levels");
    end


    % NOTE: there is a lag of 1 epoch in the timing of the radio
    if showPlots
        figure();
        yyaxis left
        plot(radio_log(1:3,:)');
        ylabel('radio on (ms)')
        ylim([0 100]) % fix scale for comparison and crop to avoid scaling caused by single outliers

        hold on;
        yyaxis right
        plot(cumsum(radio_log')/1000);

        xlabel('time (s)')
        ylabel('total radio on (s)')
        ylim([0 80]) % fix scale for comparison
        legend('FF1', 'FF2', 'FF3');

        saveFigureEps(sprintf("%s-radio", experimentNameSim));

        title("Radio on");
    end


    % take average of 3 FF for more accurate result
    total_on = (sum(sum(radio_log'))/3)/1000; %s

    samples = size(y_log, 2);

    triggers = sum(sum(((radio_log')) > 30))/3;

    trigger_ratio = triggers / samples;

    error = (y_log' - target) * 1000; % mm
    mse = sum(sum((error .^2))) / (size(error, 1) );

    ise = 0;
    iae = 0;
    itse = 0;
    iate = 0;

    for t = 1:samples
        ise = ise + sum(error(t,:) .^2);
        iae = iae + sum(abs(error(t,:)));
        itse = ise + sum(error(t,:) .^2) * t;
        iate = iae + sum(abs(error(t,:))) * t;
    end
    
    mse = ise / samples;
    mtse = itse / samples;

    %disp(experimentName);
    disp(sprintf("%s & %s & %.1f & %.2f & %.0f & %.0f  \\\\", experiments{iExperiments}(2), experiments{iExperiments}(3), total_on, trigger_ratio, mse, mtse));

        %% Load cached data
    try
        load(fileName2);
    catch
        assert(false, "Data file pstc controller does not exist");
    end
end

% % .02 9.1252e-04 (37.47)
% % 0.01 5.0639e-04 (38.31s)
% % 0.00125 4.0085e-04 (46.43)
% % per 3.9312e-04
% 
% %% Check 0,0 setting
% % convert levels to pressure value
% s1 = m2pressure(level.Data(:,1), Wis.a(3), Wis.b(3));
% s2 = m2pressure(level.Data(:,2), Wis.a(5), Wis.b(5));
% s3 = m2pressure(level.Data(:,3), Wis.a(7), Wis.b(7));
% 
% % combine sensor data, and keep only 1 sample per second
% sensor_data = [s1(1:SPS:end) s2(1:SPS:end) s3(1:SPS:end)];
% 
% % find where nothing changes changes
% sameValues = 0;
% for i = 2:samples
%     if (sensor_data(i-1, 1) == sensor_data(i, 1)) && (sensor_data(i-1, 2) == sensor_data(i, 2)) && (sensor_data(i-1, 3) == sensor_data(i, 3))
%         sameValues = sameValues + 1;
%         level.Data(i-1,1);
%         level.Data(i,1);
%         %disp('---');
%     end
% end
% 
% disp(sprintf("Same values at %d epochs, so with sigma = epsilon = 0 this would result in a trigger-ratio of: %f", sameValues, (samples - sameValues)/samples));

% minimum 