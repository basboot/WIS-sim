%% Script to compare the first and third order models using the Cantoni 
% plant parameters to show results in my thesis
% Author: Jacob Lont
% Date: 30-06-2020

%% Run the simulation and plot the first and third order responses
sim('cantonisim_distributed');
%%
t = ans.third_order.time; y = ans.third_order.signals.values;


figure(1); clf;
subplot(1,2,1)
plot(t, y,[t(1) t(end)], [1 1],'k:'); ylim([-3*10^-3 6*10^-3])
ylabel('Water level error [m]'); xlabel('time [minutes]');
legend('y5', 'y4', 'y3', 'y2', 'y1');
title('Third order plant response');
grid on;
hold on
% Plot Overshoot
% [Vmax, Imax] = max(y); plot([0 t(Imax)] ,[Vmax Vmax], 'r:')
% text(2,Vmax+0.06, ['M_p = ' num2str( round((Vmax-1)*10000)/100 ),'%'])

% Plot the settling time
% Itolarge = find(abs(y-1) > 0.02);
% if Itolarge(end) == size(y,1) 
%     % Then not settled, so dont plot
% else
%     Iset = Itolarge(end)+1;
%     Yset = y(Iset); Tset = t(Iset);
%     plot([Tset, Tset],[-0.1 Yset],'r:',[0 t(end)],[0 0],'k:')
%     text(Tset+2,Yset/2, ['t_s = ', num2str( round(100*Tset)/100 ),'s'])
% end


% hold off

% Plot the first-order response
t = ans.first_order.time; y = ans.first_order.signals.values;

% figure(1)
subplot(1,2,2)
plot(t, y,[t(1) t(end)], [1 1],'k:'); ylim([-3*10^-3 6*10^-3])
% legend('',
ylabel('Water level error [m]'); xlabel('time [minutes]');
legend('y5', 'y4', 'y3', 'y2', 'y1');
title('First order plant response');
grid on;
% hold on
% Plot Overshoot
% [Vmax, Imax] = max(y); plot([0 t(Imax)] ,[Vmax Vmax], 'r:')
% text(2,Vmax+0.06, ['M_p = ' num2str( round((Vmax-1)*10000)/100 ),'%'])

% Plot the settling time
% Itolarge = find(abs(y-1) > 0.02);
% if Itolarge(end) == size(y,1) 
%     % Then not settled, so dont plot
% else
%     Iset = Itolarge(end)+1;
%     Yset = y(Iset); Tset = t(Iset);
%     plot([Tset, Tset],[-0.1 Yset],'r:',[0 t(end)],[0 0],'k:')
%     text(Tset+2,Yset/2, ['t_s = ', num2str( round(100*Tset)/100 ),'s'])
% end


hold off