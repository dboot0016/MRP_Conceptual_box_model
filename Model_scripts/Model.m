clear all, clc
% Daan Boot, IMAU, Utrecht University, the Netherlands
% Conceptual (box) model for the Maud Rise Polynya

%% Global all necessary parameters
global rho0 rhoi Cp L sigma H h Kt Ks K alfa beta Tf FxT1 FxT2 FxS1 FxS2 NP P FT1NP FT1P FS1NP FS1P Tb Sb FT FS

%% Specify run
run='0';                                              % Run number 
                                                       % 0: normal; wn: white noise; 1a-4a, Sb, Tb, Fx, Fx1: sensitivity tests; 
n=100;                                                  % Length of simulation in years

if isequal(run,'0') || isequal(run,'wn')  % Determine which parameters are varied for sensitivity
    PM=string({'none'});
else
    PM=string({'Kt','Ks','KtKs','K','h','FF','Qoa','Qia'}); 
end

for i=1:length(PM)
    pm=PM(i);                                            % Parameter that is varied for sensitivity analysis

%% Switches for lateral fluxes (1 = on; 0 = off)
FxT1=1;                         % Lateral heat flux top layer
FxT2=1;                         % Lateral heat flux subsurface layer
FxS1=1;                         % Lateral salt flux top layer
FxS2=1;                         % Lateral salt flux subsurface layer

%% Determine which switches are on for names saved datasets
if FxT1==1 && FxT2==1 && FxS1==1 && FxS2==1         % Represents case PFB
    type='AF';                                      
elseif FxT1==1 && FxT2==1 && FxS1==1 && FxS2==0     % Represents case PFH
    type='AT';
elseif FxT1==1 && FxT2==0 && FxS1==1 && FxS2==1     % Represents case PFS
    type='AS';
elseif FxT1==0 && FxT2==1 && FxS1==0 && FxS2==1
    type='BF';
elseif FxT1==0 && FxT2==1 && FxS1==0 && FxS2==0
    type='BT';
elseif FxT1==0 && FxT2==0 && FxS1==0 && FxS2==1
    type='BS';
elseif FxT1==1 && FxT2==0 && FxS1==1 && FxS2==0
    type='CF';
elseif FxT1==1 && FxT2==0 && FxS1==0 && FxS2==0
    type='CT';
elseif FxT1==0 && FxT2==0 && FxS1==1 && FxS2==0
    type='CS';
elseif FxT1==0 && FxT2==0 && FxS1==0 && FxS2==0         % Represents case MKL and MKH
    type='DF';
end

%% Determine variables
h1=160;                         % Pycnocline depth
K1=1e-04;                       % Exchange heat coefficient between ice and 1

% Different exchange values are necessary for the different type of runs:
if isequal(type,'AF') 
    Kt1=2.82e-06;
    Ks1=2.82e-06;
end
if isequal(type,'AT') || isequal(type,'AS')
    Kt1=2.8e-06;
    Ks1=2.8e-06;   
end
if isequal(type,'CF') || isequal(type,'CT') || isequal(type,'CS') 
    Kt1=2.8e-06;
    Ks1=2.8e-06;
end
if isequal(type,'DF')
    Kt1=5e-06;          
    Ks1=1.375e-06;      % Represents MKL; change into 2e-06 for case MKH
end

%% Determine constants
AS=0.0250; fS=1/(25*365*86400); phiS=pi; DS=34.7549;            % Characteristics subsurface salt forcing
AT=0.3293; fT=1/(25*365*86400); phiT=pi; DT=0.8603;             % Characteristics subsurface heat forcing

alfa=5.82e-05; beta=8e-04;                                      % Thermal expansion, haline contraction coefficients
rhoi=900; rho0=1000; L=2.5e05; Cp=4.18e03; Tf=-1.86; sigma=30;  % Other constants

FT1NP=3.82e-9; FT1P=-3.67e-9; FS1NP=-1.35e-7; FS1P=-1.90e-7;    % Lateral fluxes top layer (determined from CESM)

dt=180; tstart=1;                                               % Time constants
H=2000;                                                         % Depth subsurface layer
FS1NP=FS1NP/22; FS1P=FS1P/22;

Tb=-0.33;                                                       % Background temperature lateral heat flux surface layer
Sb=34.4814;                                                     % Background temperature lateral salt flux surface layer
Tb1=linspace(-0.8,0.2,11);                                      % Vector for sensitivity analysis
Sb1=linspace(34.48,34.5,11);                                    % Vector for sensitivity analysis
days1=200;                                                      % Relaxation timescale
days=linspace(100,300,11);                                      % Vector for sensitivity analysis

%% Multipliers (for sensitivity analysis)
if isequal(run,'0') || isequal(run,'wn')
    x=1;
elseif isequal(run,'1a')
    x=linspace(0.5,1.5,11);
elseif isequal(run,'2a')
    x=linspace(0.8,1.2,11);
elseif isequal(run,'3a')
    x=linspace(0.1,1.9,11);
elseif isequal(run,'4a')
    x=linspace(0.9,1.1,6);
end
if isequal(pm,"Sb");
    x=1:length(Sb1);
end
if isequal(pm,"Tb");
    x=1:length(Tb1);
end
if isequal(pm,"Fx") || isequal(pm,"Fx2")
    x=1:length(days);
end

%% Determine time dependent parameters
tend=n*(86400*365);                 % End time
tsin=1:86400:(n*365*86400);         % Time vector for T2 and S2

Qoa1=load('Qoa.txt');               % Load ocean-atmosphere heat flux
Qia1=load('Qia.txt');               % Load ice-atmosphere heat flux

Qia1=circshift(Qia1,0);             % Change start time ice-atmosphere heat flux
Qoa1=circshift(Qoa1,0);             % Change start time ocean-atmosphere heat flux
Qia1=1*repmat(Qia1,1,n);            % Extent ice-atmosphere heat flux vector
Qoa1=1*repmat(Qoa1,1,n);            % Extent ocean-atmosphere heat flux vector

if isequal(run,'wn')
    NN=1:100;                       % Determine how many runs with noise
else NN=1;                          % No noise runs
end

for nn=min(NN):max(NN);
if isequal(run,'wn')                % Determine whether a white noise is run
    FF=white_noise_FF(type,n,nn);   % Function file to determine white noise
    FNP1=FF(1,:);                   % Non-polynya freshwater flux with white noise
    FP1=FF(2,:);                    % Polynya freshwater flux with white noise

else
    FNP1=load('fnp.txt');           % Load non-polynya freshwater flux
    FP1=load('fp.txt');             % Load polynya freshwater flux
    FNP1=circshift(FNP1,0);         % Change start time non-polynya freshwater flux (0 means no shift)
    FP1=circshift(FP1,0);           % Change start time polynya freshwater flux (0 means no shift)     
    FNP1=repmat(FNP1,1,n);          % Extent non-polynya freshwater flux vector
    FP1=repmat(FP1,1,n);            % Extent polynya freshwater flux vector
end

Qt=linspace(1,tend,365*n);          % Construct vector for interpolation time dependent vectors


T2=FxT2*(AT*cos(2*pi*fT*tsin+phiT))+DT*ones(1,length(tsin)); %  Lateral heat forcing subsurface layer
S2=FxS2*(AS*cos(2*pi*fS*tsin+phiS))+DS*ones(1,length(tsin)); %  Lateral salt forcing subsurface layer


%% Create empty matrices
T1=zeros(length(x),length(tsin));
d=T1; S1=T1; T2m=T1; S2m=T1; Regime=T1; 

for j=1:length(x)
%% New parameter values
if isequal(pm,"h")          %Determine which parameter is varied for sensitivity analysis
    h=x(j).*h1;
else h=h1;
end

if isequal(pm,"K") 
    K=x(j).*K1;
else K=K1;
end

if isequal(pm,"Kt")
    Kt=x(j).*Kt1;
else Kt=Kt1;
end

if isequal(pm,"Ks")
    Ks=x(j).*Ks1;
else Ks=Ks1;
end

if isequal(pm,"FF")
    FNP=x(j).*FNP1;
    FP=x(j).*FP1;
else FNP=FNP1; FP=FP1;
end

if isequal(pm,"Qoa")
    Qoa=x(j).*Qoa1;
else Qoa=Qoa1;
end

if isequal(pm,"Qia")
    Qia=x(j).*Qia1;
else Qia=Qia1;
end

if isequal(pm,"KtKs")
    Kt=x(j).*Kt1;
    Ks=x(j).*Ks1;
end

if isequal(pm,"QoaQia")
    Qoa=x(j).*Qoa1;
    Qia=x(j).*Qia1;
end

if isequal(pm,"Tb")
    Tb=Tb1(j);
else Tb=Tb;
end

if isequal(pm,"Sb")
    Sb=Sb1(j);
else Sb=Sb;
end

if isequal(pm,"Fx")
    FT=1/(days(j)*86400);
    FS=1/(days(j)*86400);
    FS1NP=1/(days(j)*86400);
    FS1P=1/(days(j)*86400);
    FT1NP=1/(days(j)*86400);
    FT1P=1/(days(j)*86400);
else FT=1/(days1*86400);
    FS=1/(days1*86400);
    FS1NP=1/(days1*86400);
    FS1P=1/(days1*86400);
    FT1NP=1/(days1*86400);
    FT1P=1/(days1*86400);
end

if isequal(pm,"Fx2")
    FT=1/(days(j)*86400);
    FS=1/(days(j)*86400);
    FS1NP=1/(days1*86400);
    FS1P=1/(days1*86400);
    FT1NP=1/(days1*86400);
    FT1P=1/(days1*86400);
else FT=1/(days1*86400);
    FS=1/(days1*86400);
    FS1NP=1/(days1*86400);
    FS1P=1/(days1*86400);
    FT1NP=1/(days1*86400);
    FT1P=1/(days1*86400);
end
tspan=tstart:dt:tend;           % Time vector

%% Initial conditions
RelT=1e-8;                       % Relative error tolerance
AbsT=1e-10;                      % Absolute error tolerance

opt=odeset('Events',@(t,y)event_4(t,y),'RelTol',RelT,'AbsTol',AbsT);
fcn=@(t,y)regime_4(t,y,Qia,Qt,FNP,FP,T2,S2);
regime=4;                                       % Initial regime

NP=1; P=0; Y=1;                                 % Start in a non-polynya period
ay=[]; t=0; T=0;                                % Renew initial conditions: values, time
ay(1,:)=[0.1 1 34.2 T2(1) S2(1)];               % Renew initial conditions: values
Regime=regime;                                  % Renew initial conditions: regime

%% Time loop for simulations
while t(end)<tspan(end);
[t,y,te,ye,ie]=ode15s(fcn, [tspan], ay(end,:),opt);

T=t(end);                                                       % Determine end time of simulation
Regime=cat(1,Regime,repmat(regime,(length(y(2:end-1,:))),1));   % Tracking of state
ay = cat(1, ay, y(2:(end-1), :));                               % Tracking of variables
ay(end,:)=y(end,:);                                             % To have right length lose one value

%% Determine next regime
if T>=tend-dt                   % If simulation stops due to end time is reached
    regime=regime;              % No regime changes
    Regime(end+1)=regime;

elseif y(end,2)<=0 && regime==4         % Regime 4 --> 2
    
    regime=2;                           % New regime
    fcn=@(t,y)regime_2(t,y,Qoa,Qt,FNP,FP,T2,S2);
    opt=odeset('Events',@(t,y)event_2(t,y),'RelTol',RelT,'AbsTol',AbsT);
    
    ay(end,2)=0;                        % Initial condition sea ice ( = 0)
    tspan=(floor(T/dt)*dt):dt:tend;     % New time span next simulation
    NP=1;
    P=0;
    
elseif y(end,2)>0 && regime==4          % Regime 4 --> 3
    
    regime=3;                           % New regime
    fcn=@(t,y)regime_3(t,y,Qia,Qt,FP,T2,S2);
    opt=odeset('Events',@(t,y)event_3(t,y),'RelTol',RelT,'AbsTol',AbsT);
    
    ay(end,1)=(h*ay(end,1)+(H-h)*ay(end,4))/(H);                    % New initial condition T: overturning
    ay(end,3)=(h*ay(end,3)+(H-h)*ay(end,5))/(H);                    % New initial condition S: overturning
    ay(end,4)=FxT2*ay(end,1)+(1-FxT2)*ay(end-1,4);                  % New initial condition T: overturning
    ay(end,5)=FxS2*ay(end,3)+(1-FxS2)*ay(end-1,5);                  % New initial condition S: overturning                           
    
    tspan=(floor(T/dt)*dt):dt:(floor(T/dt)*dt)+(1800/dt)*dt;        % New time span next simulation
    P=1; NP=0;                                                      % Switch polynya period on
     
        while Y>0 && y(end,2)>0
            regime=3;
        [t,y,te,ye,ie]=ode15s(fcn, [tspan], ay(end,:),opt);

        T=t(end);                                                   % Determine end time of simulation
        Regime=cat(1,Regime,repmat(regime,(length(y(2:end,:))),1)); % Tracking of state
        ay = cat(1, ay, y(2:(end), :));                             % Tracking of variables
        ay(end,:)=y(end,:);                                         % To have right length lose one value
    
        tspan=(floor(T/dt)*dt):dt:floor(T/dt)+(1800/dt)*dt;         % New time span next simulation
        P=1; NP=0;                                                  % Switch polynya period on
        dT=[]; dS=[];
        dT=-alfa*gradient(y(:,1));
        dS=beta*gradient(y(:,3));
        X=dT+dS;
        Y=min(X);                                                   % Variable representing static stability water column (<0 represents stable)
        end

        if y(end,2)>0
        
        regime=4;                                       % New regime
        fcn=@(t,y)regime_4(t,y,Qia,Qt,FNP,FP,T2,S2);
        opt=odeset('Events',@(t,y)event_4(t,y),'RelTol',RelT,'AbsTol',AbsT);
    
        tspan=(floor(T/dt)*dt):dt:tend;                 % New time span next simulation
         
        elseif y(end,2)<=0                                  % Regime 3 --> 1
    
        regime=1;                                       % New regime
        fcn=@(t,y)regime_1(t,y,Qoa,Qt,FP,T2,S2);
        opt=odeset('Events',@(t,y)event_1(t,y),'RelTol',RelT,'AbsTol',AbsT);                       
        
        ay(end,2)=0;
        tspan=(floor(T/dt)*dt):dt:(floor(T/dt)*dt)+(1800/dt)*dt;                 % New time span next simulation
        Y=1;
    
        while Y>0 
            regime=1;
        [t,y,te,ye,ie]=ode15s(fcn, [tspan], ay(end,:),opt);

        T=t(end);                                                       % Determine end time of simulation
        Regime=cat(1,Regime,repmat(regime,(length(y(2:end,:))),1));   % Tracking of state
        ay = cat(1, ay, y(2:(end), :));                               % Tracking of variables
        ay(end,:)=y(end,:);                                             % To have right length lose one value

        tspan=(floor(T/dt)*dt):dt:(floor(T/dt)*dt)+(1800/dt)*dt;        % New time span next simulation
        
        P=1; NP=0;                                              % Switch polynya period on
        dT=[]; dS=[];
        dT=-alfa*gradient(y(:,1));
        dS=beta*gradient(y(:,3));
        X=dT+dS;
        Y=min(X);                                               % Variable representing static stability water column (<0 means stable)
        
        end
            
        regime=2;                                       % New regime
        fcn=@(t,y)regime_2(t,y,Qoa,Qt,FNP,FP,T2,S2);
        opt=odeset('Events',@(t,y)event_2(t,y),'RelTol',RelT,'AbsTol',AbsT);
    
        tspan=(floor(T/dt)*dt):dt:tend;                 % New time span next simulation
        end

elseif y(end,1)<=Tf && regime==2                    % Regime 2 --> 4
    
    regime=4;                                       % New regime
    fcn=@(t,y)regime_4(t,y,Qia,Qt,FNP,FP,T2,S2);
    opt=odeset('Events',@(t,y)event_4(t,y),'RelTol',RelT,'AbsTol',AbsT);
    
    tspan=(floor(T/dt)*dt):dt:tend;                 % New time span next simulation
    P=0; NP=1;                                      % Switch polynya period off
    
elseif y(end,1)>Tf && regime==2                     % Regime 2 --> 1
    
    regime=1;                                       % New regime
    fcn=@(t,y)regime_1(t,y,Qoa,Qt,FP,T2,S2);
    opt=odeset('Events',@(t,y)event_1(t,y),'RelTol',RelT,'AbsTol',AbsT);

    ay(end,1)=(h*ay(end,1)+(H-h)*ay(end,4))/(H);                % New initial condition T: overturning
    ay(end,3)=(h*ay(end,3)+(H-h)*ay(end,5))/(H);                % New initial condition S: overturning
    ay(end,4)=FxT2*ay(end,1)+(1-FxT2)*ay(end-1,4);              % New initial condition T: overturning
    ay(end,5)=FxS2*ay(end,3)+(1-FxS2)*ay(end-1,5);              % New initial condition S: overturning                           
    
    tspan=(floor(T/dt)*dt):dt:(floor(T/dt)*dt)+(1800/dt)*dt;    % New time span next simulation
    P=1; NP=0;                                                  % Switch polynya period on
    Y=1; 
    
        while Y>0 
            regime=1;
        [t,y,te,ye,ie]=ode15s(fcn, [tspan], ay(end,:),opt);
        T=t(end);                                                     % Determine end time of simulation
        Regime=cat(1,Regime,repmat(regime,(length(y(2:end,:))),1));   % Tracking of state
        ay = cat(1, ay, y(2:(end), :));                               % Tracking of variables
        ay(end,:)=y(end,:);                                           % To have right length lose one value

        tspan=(floor(T/dt)*dt):dt:(floor(T/dt)*dt)+(86400/dt)*dt;     % New time span next simulation
        
        P=1; NP=0;                                                    % Switch polynya period on
        dT=[]; dS=[];
        dT=-alfa*gradient(y(:,1));
        dS=beta*gradient(y(:,3));
        X=dT+dS;
        Y=min(X);                                                     % Variable determining static stability water column (<0 means stable)
        
        end
            
        regime=2;                                       % New regime
        fcn=@(t,y)regime_2(t,y,Qoa,Qt,FNP,FP,T2,S2);
        opt=odeset('Events',@(t,y)event_2(t,y),'RelTol',RelT,'AbsTol',AbsT);
    
        tspan=(floor(T/dt)*dt):dt:tend;                 % New time span next simulation

end
end

%% Collect data
T1(j,(1:length(ay)))=ay(:,1);
d(j,(1:length(ay)))=ay(:,2);
S1(j,(1:length(ay)))=ay(:,3);
T2m(j,(1:length(ay)))=ay(:,4);
S2m(j,(1:length(ay)))=ay(:,5);
regime(j,(1:length(ay)))=(Regime(1:length(ay)))';


%% Only save daily values
ts=86400/dt;                        % Determine how many time steps in one day

for k=1:(length(T1)/ts)
    T_1(j,k)=T1(j,k*ts-ts+1);
    S_1(j,k)=S1(j,k*ts-ts+1);
    T_2(j,k)=T2m(j,k*ts-ts+1);
    S_2(j,k)=S2m(j,k*ts-ts+1);
    d_m(j,k)=d(j,k*ts-ts+1);
    
end
end

%% Save data
save(strcat(type,pm,num2str(n),run,num2str(nn)),'T_1','S_1','T_2','S_2','d_m'); % Save data

end
end
