

%% experiment
ze1 = create_iddata("20210202_step_gate1_2_s255_no_intake.csv", wis, 1, 1/128, false);
ze2 = create_iddata("20210202_step_gate1_2_s100_no_intake.csv", wis, 1, 1/128, false);
ze3 = create_iddata("20210202_step_gate1_2_s25_no_intake.csv", wis, 1, 1/128, false);

ze = merge(ze1, ze2, ze3)
%% validation

zv1 = create_iddata("20210202_step_gate1_2_s255_no_intake.csv", wis, 1, 1/128, false);
zv2 = create_iddata("20210202_step_gate1_2_s100_no_intake.csv", wis, 1, 1/128, false);
zv3 = create_iddata("20210202_step_gate1_2_s25_no_intake.csv", wis, 1, 1/128, false);

zv = merge(ze1, ze2, ze3)
%%

% create third order model with pole at origin (type 1) with a transport
% delay for identification

% pool1 1.3
% pool2 0.5 - 0.7
% pool3 1.5
init_sys = idtf(NaN(1,1),[1,NaN(1,2),0],'IODelay',0.7);

% fix last coefficient, to force tfest to keep the integrator 1/s 
init_sys.Structure.Denominator.Free = [0 1 1 0];


%%

Opt = tfestOptions('Display','on');

mtf = tfest(ze,init_sys,Opt);

%%
figure(8);
compare(zv,mtf)

