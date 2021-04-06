%% time_plot_control.m

% Show a graph to visualy check the control


% csvFile = "20210309_control_test_controller2_no_anti_windup.csv";
%% csvFile = "20210316_control_test_constant_flow5e-4.csv";
%csvFile = "20210316_control_test_controller2_no_anti_windup_min_to_sec";
 csvFile = "20210316_test_control_again";
%
%csvFile = "20210323_controller2c_offtakes_after_target";
csvFile = "20210323_controller2c";


% load data
pool_data = readmatrix(sprintf("../data/%s", csvFile));
    
figure();

plot(pool_data(:,1)/1000, pool_data(:,3)/1000);
hold on;
title("Water level pool0");

saveas(gcf,'../Latex/images/pool0', 'epsc')


figure();
plot(pool_data(:,1)/1000, pool_data(:,9)/1000);
hold on;

plot(pool_data(:,1)/1000, pool_data(:,15)/1000);
hold on;

plot(pool_data(:,1)/1000, pool_data(:,21)/1000);
hold on;
xlabel('time (s)')
ylabel('level (m)')
legend('pool1', 'pool2', 'pool3');

yline(0.25,'-','reference 1', 'LabelHorizontalAlignment', 'left', 'HandleVisibility','off');
yline(0.20,'-','reference 2', 'LabelHorizontalAlignment', 'left', 'HandleVisibility','off');
yline(0.15,'-','reference 3', 'LabelHorizontalAlignment', 'left', 'HandleVisibility','off');
title("Water levels");

saveas(gcf,'../Latex/images/pool123', 'epsc')

figure();

plot(pool_data(:,1)/1000, pool_data(:,6));
hold on;

plot(pool_data(:,1)/1000, pool_data(:,7));
hold on;

plot(pool_data(:,1)/1000, pool_data(:,6+6));
hold on;

plot(pool_data(:,1)/1000, pool_data(:,7+6));
hold on;

plot(pool_data(:,1)/1000, pool_data(:,6+12));
hold on;

plot(pool_data(:,1)/1000, pool_data(:,7+12));
hold on;
xlabel('time (s)')
ylabel('signal (x 10^{-4})')
legend('global1', 'local1', 'global2', 'local2', 'global3', 'local3');
title("Control signals");

saveas(gcf,'../Latex/images/control', 'epsc')
figure();

plot(pool_data(:,1)/1000, pool_data(:,5));
hold on;

plot(pool_data(:,1)/1000, pool_data(:,5+6));
hold on;

plot(pool_data(:,1)/1000, pool_data(:,5+12));
hold on;

xlabel('time (s)')
ylabel('servo')
legend('gate1', 'gate2', 'gate3');

title("Gates");
saveas(gcf,'../Latex/images/gates', 'epsc')
                  