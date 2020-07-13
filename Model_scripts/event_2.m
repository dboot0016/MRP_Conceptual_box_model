function [value, isterminal, direction] = event_2(t,y)
% Event function regime 2
% This function is called in Model.m when regime 2 transits to either
% regime 1 or 4

%% Global constants
global alfa beta Tf

%% Conditions
E1=((-alfa*(y(1)-y(4))+beta*(y(3)-y(5))));  % Column statically unstable (to regime 1)
E2=y(1)-Tf;                                 % Freezing temperature reached (to regime 4)

value      = [E1; E2];
isterminal = [1; 1];
direction  = [1; 0];
end