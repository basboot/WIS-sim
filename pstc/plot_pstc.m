%% 

% Time series from Simulink, with Data for 3 pools / FFs
% flow, level, servo, epoch, radio_on

target = [0.25 0.20 0.15];
showPlots = true;

useCachedData = true;

nExperiments = 1;
experiments{1} = ["20210810test0-01_1_force_trigger_python", "PSTC", "(force trigger)", "20210810test0-01_1_force_trigger_controller"];

%experiments{1} = ["20210810test0-01_1_no_trigger_python", "PSTC", "(no trigger)", "20210810test0-01_1_no_trigger_controller"];

%experiments{1} = ["20210810test0-01_1_python2", "PSTC", "(0.01 1)", "20210810test0-01_1_controller2"];

%experiments{1} = ["20210810test0-01_1_1pct_controller", "PSTC", "(0.01 1)", "20210810test0-01_1_1pct_controller"];

%experiments{1} = ["20210810test0-01_1_python_ss", "PSTC", "(0.01 1)", "20210810test0-01_1_controller_ss"];

%experiments{1} = ["20210810test0-1_2_0-1pct_python", "PSTC", "(0.1 2)", "20210810test0-1_2_0-1pct_controller"];

% possibly not the same run :-(
%experiments{1} = ["20210811test0-1_1_python", "PSTC", "(0.1 1)", "20210811test0-1_1_controller"];

experiments{1} = ["20210811test0-05_1_python", "PSTC", "(0.05 1)", "20210811test0-05_1_controller"];


% filename sim, sigma, epsilon, filename contr

for iExperiments = 1: nExperiments

    experimentNameSim = experiments{iExperiments}(1);
    experimentNameCtrl = experiments{iExperiments}(4);
    
    fileName1 = sprintf('%s.mat', experimentNameSim);
    fileName2 = sprintf('%s.mat', experimentNameCtrl);
    fileName3 = sprintf('%s.mat', "sim_0-01-1");
    %fileName3 = sprintf('%s.mat', "sim_forcetrigger");

    
    %% Load cached data
    try
        sensors = load(fileName1);
        controller = load(fileName2);
        sim = load(fileName3);
    catch
        assert(false, "Data file does not exist");
    end


    if showPlots   
        figure();
        plot(sensors.y_log');
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
        plot(sensors.radio_log');
        ylabel('radio on (ms)')
        ylim([0 100]) % fix scale for comparison and crop to avoid scaling caused by single outliers

        hold on;
        yyaxis right
        plot(cumsum(sensors.radio_log')/1000);

        xlabel('time (s)')
        ylabel('total radio on (s)')
        ylim([0 80]) % fix scale for comparison
        legend('FF1', 'FF2', 'FF3');

        saveFigureEps(sprintf("%s-radio", experimentNameSim));

        title("Radio on");
    end


    % take average of 3 FF for more accurate result
    total_on = (sum(sum(sensors.radio_log'))/3)/1000; %s

    samples = size(sensors.y_log, 2);

    triggers = sum(sum(((sensors.radio_log')) > 30))/3;

    trigger_ratio = triggers / samples;

    error = (sensors.y_log' - target) * 1000; % mm
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


    figure()
    plot(sensors.u_log');
    title('Sim - Control signal');
    
    % PSTC plots
    figure()
    plot(controller.u_log');
    ylim([0 20]);
    title('PSTC - Control signal');
    
    figure()
    plot(controller.dk_log');
    title('PSTC - Calculated sleeping times');
    
    figure()
    plot(controller.t_log');
    title('PSTC - Triggers');
    
    figure()
    plot(controller.initialized_log');
    title('PSTC - Initialisation state');
    
    figure()
    plot(controller.y_log');
    title('PSTC - Water levels');
    
    figure()
    plot(controller.radio_log');
    title('PSTC - Radio');
    
    
%     %% simulation for comparison
%     if showPlots   
%         figure();
%         plot(sim.y_log'/1000 + [0.25 0.20 0.15]);
%         hold on;
% 
%         xlabel('time (s)')
%         ylabel('level (m)')
%         legend('pool1', 'pool2', 'pool3');
% 
%         yline(0.25,'-','reference 1', 'LabelHorizontalAlignment', 'left', 'HandleVisibility','off');
%         yline(0.20,'-','reference 2', 'LabelHorizontalAlignment', 'left', 'HandleVisibility','off');
%         yline(0.15,'-','reference 3', 'LabelHorizontalAlignment', 'left', 'HandleVisibility','off');
% 
%         saveFigureEps(sprintf("%s-levels", experimentNameSim));
% 
%         title("FULL Simulation Water levels");
%         
%         figure()
%         plot(sim.u_log');
%         ylim([0 20]);
%         title('FULL SIM - Control signal');
%     end
    
  
    
end



