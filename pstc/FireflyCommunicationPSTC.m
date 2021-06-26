classdef FireflyCommunicationPSTC < handle
    %FIREFLYCOMMUNICATIONPSTC Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        port                                % comm port
        baudrate = 460800                   % comm speed                  
        device                              % serial port connection
        n_values = 13                       % num of values (doubles) per message
        log_data = []                       % log of all received data
        DEBUG = 1                           % show debug messages
        
        % PSTC
        k, xc, xptilde, X, initialized, psibar, 
        kfinal, kbar, TRIG_LEVEL, 
        np, nc, pp, mp, nw, ppt, 
        Ac, Bc, Cc, Dc, Cp, Phip, Gammap, 
        Obsbar, Vbar, V, 
        MM, Wk, QQ, Rw, Rv, wQw, cv, cvw
        
        replay_data = []
    end
    
    methods
        function obj = FireflyCommunicationPSTC(port, ...
                k, xc, xptilde, X, initialized, psibar, ...
                kfinal, kbar, TRIG_LEVEL, ...
                np, nc, pp, mp, nw, ppt, ...
                Ac, Bc, Cc, Dc, Cp, Phip, Gammap, ...
                Obsbar, Vbar, V, ...
                MM, Wk, QQ, Rw, Rv, wQw, cv, cvw)
            %FIREFLYCOMMUNICATIONPSTC Construct an instance of this class
            %   Detailed explanation goes here
            obj.port = port;
            
            obj.k = k;
            obj.xc = xc;
            obj.xptilde = xptilde;
            obj.X = X;
            obj.initialized = initialized;
            obj.psibar = psibar; 
            obj.kfinal = kfinal;
            obj.kbar = kbar;
            obj.TRIG_LEVEL = TRIG_LEVEL;
            obj.np = np;
            obj.nc = nc;
            obj.pp = pp;
            obj.mp = mp;
            obj.nw = nw;
            obj.ppt = ppt;
            obj.Ac = Ac;
            obj.Bc = Bc;
            obj.Cc = Cc;
            obj.Dc = Dc;
            obj.Cp = Cp;
            obj.Phip = Phip;
            obj.Gammap = Gammap;
            obj.Obsbar = Obsbar;
            obj.Vbar = Vbar;
            obj.V = V;
            obj.MM = MM;
            obj.Wk = Wk;
            obj.QQ = QQ;
            obj.Rw = Rw;
            obj.Rv = Rv;
            obj.wQw = wQw;
            obj.cv = cv;
            obj.cvw = cvw;
            
            
        end
        
        function connect(obj)
            % connect
            % suppress warnings about timeout, because we need this small
            % timeout to avoid blocking write
            warning ('off','all');
            obj.device = serialport(obj.port, obj.baudrate, "Timeout",0.1);
            configureTerminator(obj.device,"LF")
            
            % activate callback
            obj.activate();
        end
        
        function activate(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            configureCallback(obj.device, "byte", 4*obj.n_values, @obj.callbackMessage);
        end
        
        function deactivate(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            configureCallback(obj.device, "off");
        end
        
        function callbackMessage(obj, device, ~)
%         // epoch (1)
%         // radio on (2)
%         // send global u (3, 4, 5)
%         // send water levels (6, 7, 8, 9)
%         // send water flows (10, 11, 12)
%         // event (13)
            data = read(device, obj.n_values, "double");
            if obj.DEBUG
                disp(data(1));
            end
            obj.log_data = [obj.log_data; data];
            
            
            % handle data and send extra sleeping
            
            % multiply values by 1000 for pstc scripts
            
            % level pool 1, 2, 3 and flow 1, 2, 3
            yhat = [data(7); data(8); data(9); data(10); data(11); data(12)];
            yref = [0.25; 0.20; 0.15; 0; 0; 0];
            
            yhat = yhat - yref;
            yhat = yhat * 1000; % to mm
            
            % event detected?
            triggered = data(13);   
            % global control gate 1, 2, 3
            u = [data(3); data(4); data(5)] * 1000;
            % radio
            radio_on = data(2);
            epoch = data(1);
        
            replay = [u; yhat; triggered];
            obj.replay_data = [obj.replay_data replay];
            
           dk = obj.updatePstc(u, yhat, triggered);
            
            % force update when not initialized
            if ~obj.initialized || obj.k+2 >= obj.kfinal
                 disp('force update');
                 %obj.sendMessage("0 0"); 
            end
            
            if obj.initialized
                if (dk > 1)
                    msg = sprintf("0 %s", string(dk-1));
                    %obj.sendMessage(msg); 
                    disp("Extra sleep possible!");
                    disp(msg);
                    
                    % TODO: send sleep command and propagate the sleeping time
                end
            end
            
            fprintf("# k=%d u =(%d, %d, %d) trigger=%d, dk=%d \n", obj.k, u(1), u(2), u(3), triggered, dk);
        end
        
        function [ylog, ulog, klog] = runSimulation(obj, odeplant)
            sigma = 0.01;
            ylog = [];
            % Main simulation loop
            kk = 1;

            klog = [];
            sleeplog = [];

            %% init sim
            yref = [0.25; 0.20; 0.15; 0; 0; 0] * 1000;
            
            %% Simulation data
            x0 = 0*ones(13,1);  % Assume it's in steady state for the first moment

            % set initial water level
            % x0(1) = 0.25 * 100;
            % x0(5) = 0.2 * 100;
            % x0(9) = 0.15 * 100;
            W_MAG = -30;
            x0(end) = 0*100*W_MAG;  % This is the real disturbance
            xc0 = zeros(0,1);
            y0 = obj.Cp*x0(1:13);

            TEND = 30*60;
            obj.kfinal = 30;  % \bar{\kappa} in the paper, heartbeat = 1/2 min

            % Disturbance signal
            %omega = @(t) W_MAG*((t >= 0) & (t <= TEND/2));  % Like in NecSys

            %omega = @(t) 0.5*W_MAG*(t > 20);  % outtake after 20 secs

            %odeplant = @(t,xp,u) obj.Ap*xp + obj.Bp*u + obj.E*omega(t);


            % For reproducibility, pre-compute the noises
            rng(1907);
            V_EACH_ELEMENT = 0.1;
            h =1;
            noises = 2*V_EACH_ELEMENT*(rand(6,TEND/h + obj.kfinal + 1) - 0.5);

            %% Initialize all states need for sleep controller
            xp = x0;
            y = y0 + noises(:,1);
            yhat = y0;
            triggered = false;

            obj.k = 0;
            obj.xc = xc0;
            obj.xptilde = zeros(13,1);
            obj.X = zeros(obj.np, obj.np);
            obj.initialized = false;
            obj.psibar = zeros(obj.pp*(obj.kbar+1),1);

            sleep = 0;

            

            %% end init
            ulog = [];

            
            while kk <= TEND/h
                % Invoke controller

                % Controller
                u = obj.Cc*obj.xc + obj.Dc*yhat;
                obj.xc = obj.Ac*obj.xc + obj.Bc*yhat;

                ulog = [ulog u];

                
                % calculate sleep
                dk = obj.updatePstc(u, yhat, triggered);
%                 [dk, obj.k, obj.xc, obj.xptilde, obj.X, obj.initialized, obj.psibar] = sleepcontroller(...
%                     u, yhat, triggered, ...
%                     obj.k, obj.xc, obj.xptilde, obj.X, obj.initialized, obj.psibar, ...
%                     obj.kfinal, obj.kbar, obj.TRIG_LEVEL, ...
%                     obj.np, obj.nc, obj.pp, obj.mp, obj.nw, obj.ppt, ...
%                     obj.Ac, obj.Bc, obj.Cc, obj.Dc, obj.Cp, obj.Phip, obj.Gammap, ...
%                     obj.Obsbar, obj.Vbar, obj.V, ...
%                     obj.MM, obj.Wk, obj.QQ, obj.Rw, obj.Rv, obj.wQw, obj.cv, obj.cvw);


                if dk > 0  % There was a sleep command
                    sleep = dk;
                end
                % Run plant forward
                [tt,xpode] = ode45(@(t,x) odeplant(t, x, u), h*(kk-1:kk), xp);

                xp = xpode(end,:)';
                y = obj.Cp*xp;% + noises(:, kk+1);

                %initialized = false;

                % Trigger?  (I'm not using the sleep time here, but checking always)
                if obj.initialized
                    yh = y(1:obj.ppt);  % heights
                    yhhat = yhat(1:obj.ppt);
                    eh = yh - yhhat;
                    if any(sum(eh.^2,2) - sigma^2*sum(yh.^2,2) > obj.TRIG_LEVEL) ||...
                            obj.k+1 >= obj.kfinal
                        triggered = true;
                        yhat = y;

                        klog = [klog, obj.k];
                        sleeplog = [sleeplog, sleep];
                    else
                        triggered = false;
                    end
                else
                    yhat = y;
                end

                ylog = [ylog, yhat];

                % Iterate
                kk = kk + 1;    

            end


        end
        
        
        function ylog = runReplay(obj, odeplant)
            ylog = [];
            % Main simulation loop
            kk = 1;

            klog = [];
            sleeplog = [];

            %% sim
            yref = [0.25; 0.20; 0.15; 0; 0; 0] * 1000;
            yhat = yref * 0;
            triggered = false; % should not matter (init phase first)
            
            sigma = 0.1;

            xp = zeros(13,1);
            xp(13)=-15;
            u_log = [];
            while kk <= size(obj.replay_data, 2)
                % Invoke controller

                % Controller
                %u = obj.Cc*obj.xc + obj.Dc*yhat;
                u = obj.replay_data(1:3, kk);
                obj.xc = obj.Ac*obj.xc + obj.Bc*yhat;

%                 u_log = [u_log u];
% 
%                 % calculate sleep
%                 dk = obj.updatePstc(u, yhat, triggered);
% 
% 
%                 if dk > 0  % There was a sleep command
%                     sleep = dk;
%                 end
                % Run plant forward
                [~,xpode] = ode45(@(t,x) odeplant(t, x, u), 1*(kk-1:kk), xp);
                xp = xpode(end,:)';
                y = obj.Cp*xp;% + noises(:, kk+1);

%                 % Trigger?  (I'm not using the sleep time here, but checking always)
%                 if obj.initialized
%                     yh = y(1:3);  % heights
%                     yhhat = yhat(1:3);
%                     eh = yh - yhhat;
%                     if any(sum(eh.^2,2) - sigma^2*sum(yh.^2,2) > obj.TRIG_LEVEL) ||...
%                             obj.k+1 >= obj.kfinal
%                         triggered = true;
%                         yhat = y;
% 
%                         klog = [klog, obj.k];
%                         sleeplog = [sleeplog, sleep];
%                     else
%                         triggered = false;
%                     end
%                 else
                    yhat = y;
%                 end

                ylog = [ylog, y];

                % Iterate
                kk = kk + 1;    

            end
        end
        
        
        function dk = updatePstc(obj, u, yhat, triggered)
                [dk, obj.k, obj.xc, obj.xptilde, obj.X, obj.initialized, obj.psibar] = sleepcontroller(...
                u, yhat, triggered, ...
                obj.k, obj.xc, obj.xptilde, obj.X, obj.initialized, obj.psibar, ...
                obj.kfinal, obj.kbar, obj.TRIG_LEVEL, ...
                obj.np, obj.nc, obj.pp, obj.mp, obj.nw, obj.ppt, ...
                obj.Ac, obj.Bc, obj.Cc, obj.Dc, obj.Cp, obj.Phip, obj.Gammap, ...
                obj.Obsbar, obj.Vbar, obj.V, ...
                obj.MM, obj.Wk, obj.QQ, obj.Rw, obj.Rv, obj.wQw, obj.cv, obj.cvw);
        end
        
        function sendMessage(obj, message)
            if obj.DEBUG 
                fprintf("Send: %s\n", message);
            end
            writeline(obj.device, message)
        end
        
        function saveReplayData(obj, filename)
            replayData = obj.replay_data;
            
            save(filename, 'replayData');
        end
        
        function loadReplayData(obj, filename)
            load(filename, 'replayData');
            
            obj.replay_data = replayData;
        end        
        
    end
end

