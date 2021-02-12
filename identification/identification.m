%% identification.m

% Identify all pools

data_files = size(PoolData, 2); % total number of data files

for pool_to_identify = 1:3


    %% experiment

    clear ze

    first_iddata = true;

    for i = 1:data_files
        if PoolData(i).pool == pool_to_identify
            if PoolData(i).type == "experiment"
                if first_iddata
                    ze = createIddata(PoolData(i), Wis, false, false);
                    first_iddata = false;
                else
                    ze = merge(ze, createIddata(PoolData(i), Wis, false, false));
                end
            end
        end
    end


    %% validation

    clear zv

    first_iddata = true;

    for i = 1:data_files
        if PoolData(i).pool == pool_to_identify
            if PoolData(i).type == "validation"
                if first_iddata
                    zv = createIddata(PoolData(i), Wis, false, false);
                    first_iddata = false;
                else
                    zv = merge(zv, createIddata(PoolData(i), Wis, false, false));
                end
            end
        end
    end

    %%

    % create third order model with pole at origin (type 1) with a transport
    % delay for identification

    % pool1 1.3
    % pool2 0.5 - 0.7
    % pool3 1.5
    init_sys = idtf(NaN(1,1),[1,NaN(1,2),0],'IODelay',Wis.delays(pool_to_identify));

    % fix last coefficient, to force tfest to keep the integrator 1/s 
    init_sys.Structure.Denominator.Free = [0 1 1 0];


    %%

    Opt = tfestOptions('Display','on');
    
    clear mtf;

    mtf = tfest(ze,init_sys,Opt);

    %%
    figure();
    compare(zv,mtf)
    
    PoolModel(pool_to_identify).experiment = ze;
    PoolModel(pool_to_identify).validation = zv;
    PoolModel(pool_to_identify).tf = mtf;

end

%% Bode plots

figure();
bode(PoolModel(1).tf, PoolModel(2).tf, PoolModel(3).tf)
legend("pool1", "pool2", "pool3");

% structure of the tf is a delayed second order system combined 
% with an integrator

%                              omega_n^2                      1
%  G(s) =   e^-st * ----------------------------------  *   -----
%                   s^2 + 2 zeta omega_n s + omega_n^2        s

%% Calculate dominant wave frequency for each pool
for pool_to_identify = 1:3
    [omega_n,zeta,p] = damp(PoolModel(pool_to_identify).tf);
    % first pole is the integrator which has damping ratio -1
    assert(zeta(1) == -1, 'WARNING: something went wrong, first pole was not the integrator');

    % second and third poles are complex conjugates, just use the second
    % pole for zeta and omega)n
    PoolModel(pool_to_identify).zeta = zeta(2);
    PoolModel(pool_to_identify).omega_n = omega_n(2);
    
    % Literature survey Jacob sec. 3-3-1.
    PoolModel(pool_to_identify).phi_wave = omega_n(2) * sqrt(1-zeta(2)^2);
end
