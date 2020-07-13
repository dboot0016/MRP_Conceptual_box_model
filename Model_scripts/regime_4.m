function [dydt] = regime_4 (t,y,Qia,Qt,FNP,FP,T2,S2)
% Function file for regime 4
% This file is called on in Model.m when the model state is in regime 4

%% Global constants
global Kt K Tf rhoi L rho0 Cp h FxT1 NP FT1NP FT1P P Ks sigma FxS1 FS1NP FS1P FxT2 FxS2 H Tb Sb FT FS

%% Load time series
Qia=(interp1(Qt,Qia,t,'linear'));   % Ocean-atmosphere heat flux
FNP=(interp1(Qt,FNP,t,'linear'));   % Non-polynya freshwater flux
FP=(interp1(Qt,FP,t,'linear'));     % Polynya freshwater flux
T2=(interp1(Qt,T2,t,'linear'));     % Subsurface lateral heat forcing
S2=(interp1(Qt,S2,t,'linear'));     % Subsurface lateral salt forcing

%% ODEs
dydt=zeros(5,1);
dydt(1)=((Kt*(y(4)-y(1))-K*(y(1)-Tf))/h)+FxT1*FT1NP*(Tb-y(1))+P*FT1P*(-y(1)+Tb);                % T1
dydt(2)=(1/(rhoi*L))*(-Qia-rho0*Cp*K*(y(1)-Tf))+(NP*FNP+P*FP)/sigma;                            % Sea ice
dydt(3)=(Ks*(y(5)-y(3))-(NP*FNP+P*FP)+sigma*dydt(2))/h+FxS1*FS1NP*(Sb-y(3))+P*FS1P*(-y(3)+Sb);  % S1
dydt(4)=FxT2*(FT*(T2-y(4))+(Kt*(y(1)-y(4)))/(H-h));                                             % T2
dydt(5)=FxS2*(FS*(S2-y(5))+(Ks*(y(3)-y(5)))/(H-h));                                             % S2

end

