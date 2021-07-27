%% 

% Time series from Simulink, with Data for 3 pools / FFs
% flow, level, servo, epoch, radio_on

target = [0.25 0.20 0.15];
showPlots = false;

useCachedData = true;

nExperiments = 20;
% filename, sigma, epsilon
experiments{1} = ["continuous_disturbance_periodic", "Periodic control", ""];
experiments{2} = ["continuous_disturbance_etc_no_trigger", "ETC", "(never trigger)"];
experiments{3} = ["continuous_disturbance_etc_force_trigger", "ETC", "(always trigger)"];
experiments{4} = ["continuous_disturbance_etc_0-0-0-0", "0", "0"];
experiments{5} = ["continuous_disturbance_etc_0-1-1-flow", "0.1 (flow)", "1"];
experiments{6} = ["continuous_disturbance_etc_0-1-1", "0.1", "1"];
experiments{7} = ["continuous_disturbance_etc_0-05-1", "0.5", "1"];
experiments{8} = ["continuous_disturbance_etc_0-05-1-flow", "0.5 (flow)", "1"];
experiments{9} = ["continuous_disturbance_etc_0-01-1", "0.01", "1"];
experiments{10} = ["continuous_disturbance_etc_0-005-1", "0.005", "1"];
experiments{11} = ["continuous_disturbance_etc_0-0025-1", "0.0025", "1"];
experiments{12} = ["continuous_disturbance_etc_0-00125-1", "0.00125", "1"];
experiments{13} = ["continuous_disturbance_etc_0-01-0-1", "0.01", "0.1"];
experiments{14} = ["continuous_disturbance_etc_0-005-0-1", "0.005", "0.1"];
experiments{15} = ["continuous_disturbance_etc_0-0025-0-1", "0.0025", "0.1"]; 
experiments{16} = ["continuous_disturbance_etc_0-00125-0-1", "0.00125", "0.1"]; 
experiments{17} = ["continuous_disturbance_etc_0-01-0-01", "0.01", "0.01"];
experiments{18} = ["continuous_disturbance_etc_0-005-0-01", "0.005", "0.01"];
experiments{19} = ["continuous_disturbance_etc_0-0025-0-01", "0.0025", "0.01"];
experiments{20} = ["continuous_disturbance_etc_0-00125-0-01", "0.00125", "0.01"]; 






for iExperiments = 1: nExperiments

    experimentName = experiments{iExperiments}(1);
    
    fileName = sprintf('hil_%s.mat', experimentName);

    %% Load cached data
    try
        load(fileName);
    catch
        assert(false, "Data file does not exist");
    end


    if showPlots   
        figure();
        plot(level.Time, level.Data);
        hold on;

        xlabel('time (s)')
        ylabel('level (m)')
        legend('pool1', 'pool2', 'pool3');

        yline(0.25,'-','reference 1', 'LabelHorizontalAlignment', 'left', 'HandleVisibility','off');
        yline(0.20,'-','reference 2', 'LabelHorizontalAlignment', 'left', 'HandleVisibility','off');
        yline(0.15,'-','reference 3', 'LabelHorizontalAlignment', 'left', 'HandleVisibility','off');

        saveas(gcf,sprintf('../Latex/images/timeplot_%s',experimentName), 'epsc')
        saveFigureEps(sprintf("%s-levels", experimentName));

        title("Water levels");
    end


    % NOTE: there is a lag of 1 epoch in the timing of the radio
    if showPlots
        figure();
        yyaxis left
        plot(radio_on.Time(1:SPS:end), double(radio_on.Data(1:SPS:end, :)));
        ylabel('radio on (ms)')
        ylim([0 100]) % fix scale for comparison and crop to avoid scaling caused by single outliers

        hold on;
        yyaxis right
        plot(radio_on.Time(1:SPS:end), cumsum(double(radio_on.Data(1:SPS:end, :)))/1000);

        xlabel('time (s)')
        ylabel('total radio on (s)')
        ylim([0 80]) % fix scale for comparison
        legend('FF1', 'FF2', 'FF3');

        saveas(gcf,sprintf('../Latex/images/radio_%s',experimentName), 'epsc')
        saveFigureEps(sprintf("%s-radio", experimentName));

        title("Radio on");
    end


    % take average of 3 FF for more accurate result
    total_on = (sum(sum(double(radio_on.Data(1:SPS:end, :))))/3)/1000; %s

    samples = size(radio_on.Data(1:SPS:end, :), 1);

    triggers = sum(sum(((double(radio_on.Data(1:SPS:end, :)))) > 30))/3;

    trigger_ratio = triggers / samples;

    error = (level.Data - target);
    error = error(1:SPS:end, :); % only keep error at sampling times
    mse = sum(sum((error .^2))) / (size(error, 1) );

    ise = 0;
    iae = 0;
    itse = 0;
    iate = 0;

    for t = 1:samples
        ise = ise + sum(error(i,:) .^2);
        iae = iae + sum(abs(error(i,:)));
        itse = ise + sum(error(i,:) .^2) * t;
        iate = iae + sum(abs(error(i,:))) * t;
    end

    %disp(experimentName);
    disp(sprintf("%s & %s & %f & %f & %f & %f & %f & %f \\\\", experiments{iExperiments}(2), experiments{iExperiments}(3), total_on, trigger_ratio, ise, iae, itse, iate));

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