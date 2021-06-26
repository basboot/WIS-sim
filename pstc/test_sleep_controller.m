% Test sleep_controller.m

clear all;
init_plant;
init_etc;
init_pstc;
init_simulation;

ylog = [];
% Main simulation loop
kk = 1;

klog = [];
sleeplog = [];

while kk <= TEND/h
    % Invoke controller
    [u, dk, k, xc, xptilde, X, initialized, psibar] = sleepcontroller(...
        yhat, triggered, ...
        k, xc, xptilde, X, initialized, psibar, ...
        kfinal, kbar, TRIG_LEVEL, ...
        np, nc, pp, mp, nw, ppt, ...
        Ac, Bc, Cc, Dc, Cp, Phip, Gammap, ...
        Obsbar, Vbar, V, ...
        MM, Wk, QQ, Rw, Rv, wQw, cv, cvw);
    
    if dk > 0  % There was a sleep command
        sleep = dk;
    end
    % Run plant forward
    [tt,xpode] = ode45(@(t,x) odeplant(t, x, u), h*(kk-1:kk), xp);
    xp = xpode(end,:)';
    y = Cp*xp + noises(:, kk+1);
    
    % Trigger?  (I'm not using the sleep time here, but checking always)
    if initialized
        yh = y(1:ppt);  % heights
        yhhat = yhat(1:ppt);
        eh = yh - yhhat;
        if any(sum(eh.^2,2) - sigma^2*sum(yh.^2,2) > TRIG_LEVEL) ||...
                k+1 >= kfinal
            triggered = true;
            yhat = y;
            
            klog = [klog, k];
            sleeplog = [sleeplog, sleep];
        else
            triggered = false;
        end
    else
        yhat = y;
    end
    
    ylog = [ylog, y];
            
    % Iterate
    kk = kk + 1;    
    
end