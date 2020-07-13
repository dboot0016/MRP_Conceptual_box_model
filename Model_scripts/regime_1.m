function [dydt] = regime_1 (t,y,Qoa,Qt,FP,T2,S2)
% Function file for regime 1
% This file is called on in Model.m when the model state is in regime 1

%% Global constants
global rho0 Cp h FxT1 FT1P FxS1 FS1P FxT2 FxS2 H Tb Sb FT FS

%% Load time series
Qoa=(interp1(Qt,Qoa,t,'linear'));   % Ocean-atmosphere heat flux
FP=(interp1(Qt,FP,t,'linear'));     % Polynya freshwater flux
T2=(interp1(Qt,T2,t,'linear'));     % Subsurface lateral heat forcing
S2=(interp1(Qt,S2,t,'linear'));     % Subsurface lateral salt forcing

%% ODEs
dydt=zeros(5,1);
dydt(1)=(Qoa./(rho0*Cp)+FxT1*FT1P*(-y(1)+Tb)*h)/H+FxT2*FT*(T2-y(1))*(H-h)/H;    % T1
dydt(2)=0;                                                                      % Sea ice
dydt(3)=(-FP+FxS1*FS1P*(-y(3)+Sb)*h)/H+FxS2*FS*(S2-y(3))*(H-h)/H;               % S1
dydt(4)=FxT2*dydt(1);                                                           % T2
dydt(5)=FxS2*dydt(3);                                                           % S2

end

