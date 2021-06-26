%% 

% Time series from Simulink, with Data for 3 pools / FFs
% flow, level, servo, epoch, radio_on

target = [0.25 0.20 0.15];

useCachedData = true;
%experimentName = "continuous_disturbance_etc_0-00125";
%experimentName = "continuous_disturbance_etc_0-0025";
%experimentName = "continuous_disturbance_etc_0-005";
%experimentName = "continuous_disturbance_etc_0-01"; 
%experimentName = "continuous_disturbance_etc_0-02"; 

experimentName = "continuous_disturbance_etc_new_correct_0-01";
fileName = sprintf('hil_%s.mat', experimentName);

%% Load cached, or save new data
if useCachedData
    try
        load(fileName);
    catch
        assert(false, "Data file does not exist");
    end
else
    save(fileName,'flow','level','servo','epoch',...
        'radio_on', 'SPS');    
end

   
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



% NOTE: there is a lag of 1 epoch in the timing of the radio
figure();
yyaxis left
plot(radio_on.Time(1:SPS:end), double(radio_on.Data(1:SPS:end, :)));
ylabel('radio on (ms)')

hold on;
yyaxis right
plot(radio_on.Time(1:SPS:end), cumsum(double(radio_on.Data(1:SPS:end, :)))/1000);

sum(double(radio_on.Data(1:SPS:end, :)))

xlabel('time (s)')
ylabel('total radio on (s)')
legend('FF1', 'FF2', 'FF3');

saveas(gcf,sprintf('../Latex/images/radio_%s',experimentName), 'epsc')
saveFigureEps(sprintf("%s-radio", experimentName));

title("Radio on");

error = level.Data - target;
mse = sum(sum((error .^2))) / size(error, 1) * size(error, 2)

% .02 9.1252e-04 (37.47)
% 0.01 5.0639e-04 (38.31s)
% 0.00125 4.0085e-04 (46.43)
% per 3.9312e-04