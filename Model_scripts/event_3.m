function [value, isterminal, direction] = event_3(t,y)
% Event function regime 3
% This function is called in Model.m.
% It is only necessary for the ODE15s solver, what's in this script
% exactly, is not important.
% The regime transitions in Model.m for regimes 1 and 3 are entered manually.

%% Global constants

%% Conditions
E1=1;
E2=1;

value      = [E1; E2];
isterminal = [1; 1];
direction  = [1; 0];
end