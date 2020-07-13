function [value, isterminal, direction] = event_4(t,y)
% Event function regime 4
% This function is called in Model.m when regime 4 transits to either
% regime 2 or 3

%% Global constants
global alfa beta

%% Conditions
E1=-alfa*(y(1)-y(4))+beta*(y(3)-y(5));  % Column statically unstable (to regime 3)
E2=y(2);                                % Sea ice = 0 (to regime 2)

value      = [E1; E2];
isterminal = [1; 1];
direction  = [1; 0];
end