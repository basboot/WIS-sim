function [flow] = gate2flow(servo, h1, h2)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% flow in m^3/sec

K = 6.5000e-05 *60; % controller asks for flow per minute so gate constant must be multiplied
delta_level = h1 - h2;

% calculate flow with actual gate setting
flow = (servo + 0) * (K )  * sqrt(abs(delta_level)) * sign(delta_level); 

