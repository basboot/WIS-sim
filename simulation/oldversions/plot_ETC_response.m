%% Script to plot ETC controlled simulink model responses 
% using plant parameters to show results in my thesis
% Author: Jacob Lont
% Date: 04-08-2020

%% Run the simulation and plot the first and third order responses
sim('cantoni_combined_matrices_ETC');
out = ans;
%% Plot the signals
y = out.outputs_cont;       % Continuously measured
u = out.inputs_cont;        % Continuously measured
y_etc = out.outputs_etc;    % ETC signals
u_etc = out.inputs_etc;     % ETC signals
triggers = out.triggers;    % Triggering instances
tmax = max(out.tout);       % Equal to the number of samples using h=1
disturbance_yref = out.dist_yref;    % The reference set-point disturbances

y_lim_min = -0.12;
y_lim_max = 0.2;

u_lim_min = -0.12;
u_lim_max = 0.1;

% output disturbances
d_lim_min = -0.05;
d_lim_max = 0.15;

% flow disturbances
d_lim_min_u = -0.05;
d_lim_max_u = 0.15;

g = figure(1); clf;
g.Resize = 'off';
set(g, 'Position',  [100, 100, 800, 800])
numplots = 4;

% Plot the plant outputs signals (ETC)
subplot(numplots,1,1)
plot(y_etc); %plot the timeseries
ylim([y_lim_min, y_lim_max]);
ylabel('Water level [m]'); xlabel('time [minutes]');
legend('y1', 'y2', 'y3', 'y4', 'y5');
title(strcat('Plant outputs using ETC w.r.t. initial setpoints for $\sigma$ = ', num2str(sigma)), 'interpreter','latex');
grid on;
hold on

% Plot the control signals (ETC)
subplot(numplots,1,2)
plot(u_etc); %plot the timeseries
ylim([u_lim_min, u_lim_max]);
ylabel('Control signal [m^3/min]'); xlabel('time [minutes]');
legend('u1', 'u2', 'u3', 'u4', 'u5');
title('Control signals in ETC', 'interpreter', 'latex');
grid on;

% Plot the set-point reference disturbance profile
subplot(numplots,1,3)
plot(disturbance_yref); %plot the timeseries
ylim([d_lim_min, d_lim_max]);
ylabel('Set-point reference [m]'); xlabel('time [minutes]');
legend('r1', 'r2', 'r3', 'r4', 'r5');
title('Simulated disturbance profile', 'interpreter', 'latex');
grid on;


% Plot the ETC triggering instances
subplot(numplots,1,4)
plot(triggers,'s-b'); %plot the timeseries
ylabel('Trigger'); xlabel('time [minutes]');
yticks([0 1])
yticklabels({'False','True'})
triggercolumn = out.triggers.data;
title(strcat(num2str(sum(triggercolumn)), {' '}, 'triggers of', {' '}, num2str(tmax), {' '}, 'sample instances'), 'interpreter', 'latex');
grid on;


hold off

%% Save the figure using a command
% saveas(g, 'figures\ETCplot_udist_0p10.eps','epsc'); close;



