clear all,clc

% Daan Boot, IMAU, Utrecht University
% Script used to analyse salinity of CESM run
% Part 1: Layer 1 (0-160m) 
% Part 2: Layer 2 (160-2000m)
% Part 3: Lateral flux layer 1 (0-160m)
% End product(s): A fit for subsurface salt forcing;
%                 Values for the lateral salt flux in layer 1.  
% CESM available upon request

%% Load in CESM
S=ncread('SALT_Polynya_depth.nc','SALT');   % CESM data
S=S(:,602:1813);                            % CESM data years 150-250
l1=3; l2=12;                                % Indices for first layer (depth)
l3=14; l4=22;                               % Indices for second layer (depth)

S1=S(l1:l2,:);                              % Temperature layer 1
S2=S(l3:l4,:);                              % Temperature layer 2

l=ncread('Layer_depth.nc','layer');         % Size per layer [m]
d=ncread('Layer_depth.nc','depth');         % Mean depth at cell
D=-d(1:41);                                 % Mean depth converted for plotting

t=linspace(150,250,length(S1));             % Time axis for plotting
t11=121; t12=396;                           % Time indices period 1 (after polynya 1)
t21=397; t22=708;                           % Time indices period 2 (after polynya 2)
t31=709; t32=1032;                          % Time indices period 3 (after polynya 3)
%% Part 1: Top layer 1

%% Contour plots : Temperature vs depth for the three periods
figure(1)
contourf(t(t11:t12),D(l1:l2),S1(:,t11:t12))
grid on
box on
xlabel('Time [CESM model years]','FontSize',12)
ylabel('Depth [m]','FontSize',12)
title('S_1 for first period (160-182)','FontSize',12)
colorbar

figure(2)
contourf(t(t21:t22),D(l1:l2),S1(:,t21:t22))
grid on
box on
xlabel('Time [CESM model years]','FontSize',12)
ylabel('Depth [m]','FontSize',12)
title('S_1 for second period (183-208)','FontSize',12)
colorbar

figure(3)
contourf(t(t31:t32),D(l1:l2),S1(:,t31:t32))
grid on
box on
xlabel('Time [CESM model years]','FontSize',12)
ylabel('Depth [m]','FontSize',12)
title('S_1 for third period (209-235)','FontSize',12)
colorbar

%% Line plots: Mean temperature vs time
S1m=sum((S1(:,:)).*l(l1:l2))/sum(l(l1:l2));

figure(4)
plot(t(t11:t32),S1m(t11:t32),'LineWidth',1.5)
xlabel('Time [CESM model years]','FontSize',12)
ylabel('S_1 [^{\circ} C]','FontSize',12)
title('Mean temperature top layer','FontSize',12)
xlim([t(t11) t(t32)])
grid on
box on

%% Part 2: Subsurface layer (160-2000)

%% Contour plots : Temperature vs depth for the three periods
figure(5)
contourf(t(t11:t12),D(l3:l4),S2(:,t11:t12))
grid on
box on
xlabel('Time [CESM model years]','FontSize',12)
ylabel('Depth [m]','FontSize',12)
title('S_2 for first period (160-182)','FontSize',12)
colorbar

figure(6)
contourf(t(t21:t22),D(l3:l4),S2(:,t21:t22))
grid on
box on
xlabel('Time [CESM model years]','FontSize',12)
ylabel('Depth [m]','FontSize',12)
title('S_2 for second period (183-208)','FontSize',12)
colorbar

figure(7)
contourf(t(t31:t32),D(l3:l4),S2(:,t31:t32))
grid on
box on
xlabel('Time [CESM model years]','FontSize',12)
ylabel('Depth [m]','FontSize',12)
title('S_2 for third period (209-235)','FontSize',12)
colorbar

%% Line plots: Mean temperature vs time
S2m=sum((S2(:,:)).*l(l3:l4))/sum(l(l3:l4));

figure(8)
plot(t(t11:t32),S2m(t11:t32),'LineWidth',1.5)
xlabel('Time [CESM model years]','FontSize',12)
ylabel('S_2 [g/kg]','FontSize',12)
title('Mean salinity subsurface layer','FontSize',12)
xlim([t(t11) t(t32)])
grid on
box on

%% Fit subsurface forcing to data
fitS=S2m(t21-12:t22-12);
a=0.5*(max(fitS)-(min(fitS)+0.01))
b=(1/(25*86400*365));
c=pi;
d=0.5*(max(fitS)+(min(fitS)+0.01))

tsin=1:86400:25*365*86400;
SF=a*cos(b*2*pi*tsin+c)+d;
fitS=interp1(linspace(1,length(tsin),length(fitS)),fitS,1:length(tsin)); % Eventual fit for subsurface salt flux

xx=260;
yy=200;
FS=8;

figure(9)
hold on
plot(tsin/(86400*365),fitS,'LineWidth',1.5)
plot(tsin/(86400*365),SF,'LineWidth',1.5)
xlabel('Time [years]','FontSize',FS)
ylabel('Salinity [g/kg]','FontSize',FS)
title('Fitted S_{b2}','FontSize',FS)
legend({'CESM data','Fitted S_{b2}'},'Location','south','FontSize',FS)
grid on
box on
ax = gca;
ax.FontSize = FS; 
xlim([0 25])
ylim([34.715 34.785])

x0=25;
y0=150;
width=xx;
height=yy;
set(gcf,'position',[x0,y0,width,height])

%print(gcf,'figure11(b).png','-dpng','-r1200');

%% Part 3: Lateral flux layer 1 (0-160m)

%% Load in variables CESM
usaltE=ncread('Salt_flux_east.nc','USALT');     % Salt flux through eastern border
usaltW=ncread('Salt_flux_west.nc','USALT');     % Salt flux through western border
vsaltN=ncread('Salt_flux_north.nc','VSALT');    % Salt flux through northern border
vsaltS=ncread('Salt_flux_south.nc','VSALT');    % Salt flux through southern border

%% Transform data
usaltE=usaltE(:,541:1872);                      % Select time period
usaltW=usaltW(:,541:1872);                      % Select time period
vsaltN=vsaltN(:,541:1872);                      % Select time period
vsaltS=vsaltS(:,541:1872);                      % Select time period

l11=l1; l21=l2;                                 % To be able to vary depth
l31=l3; l41=l4;                                 % To be able to vary depth

S1E=sum(usaltE(l11:l21,:));                     % Sum over layers [layer 1]
S1W=sum(usaltW(l11:l21,:));                     % Sum over layers [layer 1]
S1N=sum(vsaltN(l11:l21,:));                     % Sum over layers [layer 1]
S1S=sum(vsaltS(l11:l21,:));                     % Sum over layers [layer 1]

S1F=-S1E+S1W-S1N+S1S;                           % Determine total flux layer 1
S1Fm=movmean(S1F,60);                           % Take moving average over 60 months
S1Fm=S1Fm(61:end-60);                           % Select right time period (150-250) 

%% Plot flux in period 2
figure(10)
plot(t(t21:t22),S1Fm(t21:t22),'LineWidth',1.5)
grid on
box on
xlabel('Time [CESM model years]','FontSize',12)
ylabel('Salt flux [psu Sv]','FontSize',12)
title('Lateral salt flux layer 1 for period 2','FontSize',12)

%% Determine values flux
t21a=t22-48;                                    % End time non-polynya period
S1NP=mean(S1Fm(t21:t21a));                      % Take mean of non-polynya period
S1P=mean(S1Fm(t21a+1:t22));                     % Take mean of polynya period

rho=1026;                                       % Density of sea water
S1NP=S1NP*rho*1e06*1e-03;                       % Convert to kg/s
S1P=S1P*rho*1e06*1e-03;                         % Convert to kg/s

A=1.44e11;                                      % Area of polynya region
h=160;                                          % Depth of layer 1
salt=34.5e-03;                                  % Salinity of layer 1
TS=A*h*rho*salt;                                % Total salt present in layer 1 [kg]

S_altNP=S1NP/TS*salt*1e03;                      % Determine salt flux using ratio and salinity (non-polynya)
S_altP=S1P/TS*salt*1e03;                        % Determine salt flux using ratio and salinity (polynya)

S1NP=S1NP*1e03/(rho*A*h)                        % Determine salt flux by converting kg/s to g/kg/s
S1P=S1P*1e03/(rho*A*h)                          % Determine salt flux by converting kg/s to g/kg/s

%% Background salinity CESM
Polynya=[337:396 661:708 973:1032];                     % Determine polynya periods
N=1:1212;                                               % Entire time vector
N(Polynya)=[]; 

S1NPm=S1m(N);

mean(S1NPm)                                             % Fitted restoring value for lateral salt flux surface layer
