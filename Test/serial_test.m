clear all;

port = "/dev/cu.SLAB_USBtoUART10";
baudrate = 460800;

device = serialport(port, baudrate);


% s.BytesAvailableFcnMode = 'terminator';
% s.BytesAvailableFcn = @instrcallback;

%         // epoch (1)
%         // radio on (2)
%         // send global u (3, 4, 5)
%         // send water levels (6, 7, 8, 9)
%         // send water flows (10, 11, 12)

flush(device)
configureTerminator(device,"LF")

DEBUG = 0;


if DEBUG
    configureCallback(device, "terminator" ,@callbackLogging)
else
    configureCallback(device, "byte", 4*12, @callbackDouble)
end

% while 1
%     disp('next')
%     readline(device)
% end
% 
% fopen(s)

function callbackLogging(device, event)
    disp('event');
    data = readline(device);
    disp(data);
end

function callbackDouble(device, event)
    %disp('event');
    %data = readline(device);
    data = read(device,12,"double");
    %disp(data);
    %C = strsplit(data);
%     for i = 1:size(data, 1)
%         disp(C(i));
%     end
    %values = double(char(data));
    % skip all debug info
    %if values(1) ~= 91
    %    data
    %end
    
    data(3) = data(3) / 10000;
    disp(data);
%            data = fgetl(app.serial_handler);
           % ...
end

%writeline(device,sprintf("%s %s", string(247), string(0)))

% configureCallback(device,"off")