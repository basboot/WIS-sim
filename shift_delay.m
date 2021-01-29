function [input,output] = shift_delay(measured_input,measured_output,delay)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[M, ~] = size(measured_input);

input = measured_input(1:end-delay);
output = measured_output(delay+1:end);

end

