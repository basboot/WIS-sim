clear all;

port = "/dev/cu.SLAB_USBtoUART3";
baudrate = 460800;

device = serialport(port, baudrate);

flush(device)
readline(device)
