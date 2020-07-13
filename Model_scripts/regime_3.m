function [dydt] = regime_3 (t,y,Qia,Qt,FP,T2,S2)
% Function file for regime 3
% This file is called on in Model.m when the model state is in regime 3

%% Global constants
global K Tf rhoi L rho0 Cp h FxT1 FT1P sigma FxS1 FS1P FxT2 FxS2 H Tb Sb FT FS

%% Load time series
Qia=(interp1(Qt,Qia,t,'linear'));   % Ocean-atmosphere heat flux
FP=(interp1(Qt,FP,t,'linear'));     % Polynya freshwater flux
T2=(interp1(Qt,T2,t,'linear'));     % Subsurface lateral heat forcing
S2=(interp1(Qt,S2,t,'linear'));     % Subsurface lateral salt forcing

%% ODEs
dydt=zeros(5,1);
dydt(1)=(K*(y(1)-Tf)+(FxT1*FT1P*h*(-y(1)+Tb)))/H+FxT2*FT*(T2-y(1))*(H-h)/H;     % T1
dydt(2)=(1/(rhoi*L))*(-Qia-rho0*Cp*K*(y(1)-Tf))+FP/sigma;                       % Sea ice
dydt(3)=(-FP+sigma*dydt(2)+FxS1*FS1P*h*(-y(3)+Sb))/H+FxS2*FS*(S2-y(3))*(H-h)/H; % S1
dydt(4)=FxT2*dydt(1);                                                           % T2
dydt(5)=FxS2*dydt(3);                                                           % S2

end

