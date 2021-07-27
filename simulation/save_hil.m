%experimentName = "continuous_disturbance_etc_0-01-1";
% experimentName = "continuous_disturbance_etc_0-005-1";
% experimentName = "continuous_disturbance_etc_0-0025-1";
% experimentName = "continuous_disturbance_etc_0-00125-1";

%experimentName = "continuous_disturbance_etc_0-01-0-1";
%experimentName = "continuous_disturbance_etc_0-005-0-1";
%experimentName = "continuous_disturbance_etc_0-0025-0-1"; 
%experimentName = "continuous_disturbance_etc_0-00125-0-1"; 

%experimentName = "continuous_disturbance_etc_0-0025-0-01";
%experimentName = "continuous_disturbance_etc_0-005-0-01";
%experimentName = "continuous_disturbance_etc_0-01-0-01";
%experimentName = "continuous_disturbance_etc_0-00125-0-01"; 

%experimentName = "continuous_disturbance_etc_0-0-0-0";

%experimentName = "continuous_disturbance_etc_force_trigger";

%experimentName = "continuous_disturbance_etc_no_trigger";

%experimentName = "continuous_disturbance_periodic";

%experimentName = "continuous_disturbance_etc_0-05-1";

%experimentName = "continuous_disturbance_etc_0-05-1-flow";

%experimentName = "continuous_disturbance_etc_0-1-1-flow";

experimentName = "continuous_disturbance_etc_0-1-1";


fileName = sprintf('hil_%s.mat', experimentName);

save(fileName,'flow','level','servo','epoch',...
        'radio_on', 'SPS');    
