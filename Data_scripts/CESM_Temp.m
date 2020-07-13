clear all,clc

% Daan Boot, IMAU, Utrecht University
% Script used to analyse temperature of CESM run
% Part 1: Layer 1 (0-160m) 
% Part 2: Layer 2 (160-2000m)
% Part 3: Lateral flux in layer 1 (0-160m)
% End product(s): A fit for subsurface heat forcing
%                 Values for the lateral heat flux in layer 1
% CESM data available upon request

%% Load in CESM
T=ncread('TEMP_Polynya_depth.nc','TEMP');   % CESM data
T=T(:,602:1813);                            % CESM data years 150-250
l1=3; l2=12;                                % Indices for first layer (depth)
l3=14; l4=22;                               % Indices for second layer (depth)

T1=T(l1:l2,:);                              % Temperature layer 1
T2=T(l3:l4,:);                              % Temperature layer 2

l=ncread('Layer_depth.nc','layer');         % Size per layer [m]
d=ncread('Layer_depth.nc','depth');         % Mean depth at cell
D=-d(1:41);                                 % Mean depth converted for plotting

t=linspace(150,250,length(T1));             % Time axis for plotting
t11=121; t12=396;                           % Time indices period 1 (after polynya 1)
t21=397; t22=708;                           % Time indices period 2 (after polynya 2)
t31=709; t32=1032;                          % Time indices period 3 (after polynya 3)
%% Part 1: Top layer 1

%% Contour plots : Temperature vs depth for the three periods
figure(1)
contourf(t(t11:t12),D(l1:l2),T1(:,t11:t12))
grid on
box on
xlabel('Time [CESM model years]','FontSize',12)
ylabel('Depth [m]','FontSize',12)
title('T_1 for first period (160-182)','FontSize',12)
colorbar

figure(2)
contourf(t(t21:t22),D(l1:l2),T1(:,t21:t22))
grid on
box on
xlabel('Time [CESM model years]','FontSize',12)
ylabel('Depth [m]','FontSize',12)
title('T_1 for second period (183-208)','FontSize',12)
colorbar

figure(3)
contourf(t(t31:t32),D(l1:l2),T1(:,t31:t32))
grid on
box on
xlabel('Time [CESM model years]','FontSize',12)
ylabel('Depth [m]','FontSize',12)
title('T_1 for third period (209-235)','FontSize',12)
colorbar

%% Line plots: Mean temperature vs time
T1m=sum((T1(:,:)).*l(l1:l2))/sum(l(l1:l2));

figure(4)
plot(t(t11:t32),T1m(t11:t32),'LineWidth',1.5)
xlabel('Time [CESM model years]','FontSize',12)
ylabel('T_1 [^{\circ} C]','FontSize',12)
title('Mean temperature top layer','FontSize',12)
xlim([t(t11) t(t32)])
grid on
box on

%% Part 2: Subsurface layer (160-2000)

%% Contour plots : Temperature vs depth for the three periods
figure(5)
contourf(t(t11:t12),D(l3:l4),T2(:,t11:t12))
grid on
box on
xlabel('Time [CESM model years]','FontSize',12)
ylabel('Depth [m]','FontSize',12)
title('T_2 for first period (160-182)','FontSize',12)
colorbar

figure(6)
contourf(t(t21:t22),D(l3:l4),T2(:,t21:t22))
grid on
box on
xlabel('Time [CESM model years]','FontSize',12)
ylabel('Depth [m]','FontSize',12)
title('T_2 for second period (183-208)','FontSize',12)
colorbar

figure(7)
contourf(t(t31:t32),D(l3:l4),T2(:,t31:t32))
grid on
box on
xlabel('Time [CESM model years]','FontSize',12)
ylabel('Depth [m]','FontSize',12)
title('T_2 for third period (209-235)','FontSize',12)
colorbar

%% Line plots: Mean temperature vs time
T2m=sum((T2(:,:)).*l(l3:l4))/sum(l(l3:l4));

figure(8)
plot(t(t11:t32),T2m(t11:t32),'LineWidth',1.5)
xlabel('Time [CESM model years]','FontSize',12)
ylabel('T_2 [^{\circ} C]','FontSize',12)
title('Mean temperature subsurface layer','FontSize',12)
xlim([t(t11) t(t32)])
grid on
box on

%% Fit subsurface forcing to data
fitT=T2m(t21-12:t22-12);                    % Select period to be fitted
a=0.5*(max(fitT)-(min(fitT)+0.1))           % Amplitude fit
b=(1/(25*86400*365));                       % Frequency fit
c=pi;                                       % Phase fit
d=0.5*(max(fitT)+(min(fitT)+0.1))           % Offset fit

tsin=1:86400:25*365*86400;                  % Time vector for fit
SF=a*cos(b*2*pi*tsin+c)+d;                  % Subsurface forcing (temperature)
fitT=interp1(linspace(1,length(tsin),length(fitT)),fitT,1:length(tsin));    % Interpolate CESM to be able to plot 

xx=260;
yy=200;
FS=8;

figure(9)
hold on
plot(tsin/(86400*365),fitT,'LineWidth',1.5)
plot(tsin/(86400*365),SF,'LineWidth',1.5)
xlabel('Time [years]','FontSize',FS)
ylabel('Temperature [^{\circ} C]','FontSize',FS)
title('Fitted T_{b2}','FontSize',FS)
legend({'CESM data','Fitted T_{b2}'},'Location','south','FontSize',FS)
grid on
box on
xlim([0 25])
ax = gca;
ax.FontSize = FS; 

x0=25;
y0=150;
width=xx;
height=yy;
set(gcf,'position',[x0,y0,width,height])

print(gcf,'figure11(a).png','-dpng','-r1200');

%% Part 3: Lateral flux layer 1 (0-160m)

%% Load in CESM variables
uheatE=ncread('Heat_flux_east.nc','UHEAT');
uheatW=ncread('Heat_flux_west.nc','UHEAT');
vheatN=ncread('Heat_flux_north.nc','VHEAT');
vheatS=ncread('Heat_flux_south.nc','VHEAT');

%% Transform data
uheatE=uheatE(:,541:1872);                      % Select right time period
uheatW=uheatW(:,541:1872);                      % Select right time period
vheatN=vheatN(:,541:1872);                      % Select right time period
vheatS=vheatS(:,541:1872);                      % Select right time period

l11=l1; l21=l2;                                 % To be able to vary depth layer 1
l31=l3; l41=l4;                                 % To be able to vary depth layer 2

T1E=sum(uheatE(l11:l21,:));                     % Sum over depth layer 1
T1W=sum(uheatW(l11:l21,:));                     % Sum over depth layer 1
T1N=sum(vheatN(l11:l21,:));                     % Sum over depth layer 1
T1S=sum(vheatS(l11:l21,:));                     % Sum over depth layer 1

T1F=-T1E+T1W-T1N+T1S;                           % Determine total flux
T1Fm=movmean(T1F,60);                           % Take moving average over 60 months
T1Fm=T1Fm(61:end-60);                           % Select years 150-250;

%% Plot flux in period 2
figure(10)
plot(t(t21:t22),T1Fm(t21:t22),'LineWidth',1.5)
grid on
box on
xlabel('Time [CESM model years]','FontSize',12)
ylabel('Heat flux [W]','FontSize',12)
title('Lateral heat flux layer 1 for period 2','FontSize',12)

%% Determine values
t21a=t22-48;                                    % End time non-polynya period
T1NP=mean(T1Fm(t21:t21a));                      % Take mean of non-polynya period
T1P=mean(T1Fm(t21a+1:t22));                     % Take mean of polynya period

rho=1000;                                       % Density water
A=1.44e11;                                      % Area polynya region
h=160;                                          % Depth layer 1
Cp=4.18e03;                                     % Specific heat of seawater

T1NP=T1NP/(rho*A*h*Cp)                          % Transform W to K/s
T1P=T1P/(rho*A*h*Cp)                            % Transform W to K/s

%% Background temperature CESM
Polynya=[337:396 661:708 973:1032];                     % Determine polynya periods
N=1:1212;                                               % Entire time vector
N(Polynya)=[]; 

T1NPm=T1m(N);

mean(T1NPm)                                             % Fit for the lateral heat flux in the surface layer

%%
points=[t(661) t(661) t(708) t(708)];
points2=[0.3 1.2 1.2 0.3];
figure(1)
hold on
plot(t(409:708),T2m(409:708),'LineWidth',1.5)
patch(points,points2,'black','FaceAlpha',0.3)

grid on
box on
xlim([t(408) t(708)])
ylim([0.3 1.2])
xlabel('Time [CESM model years]','FontSize',12)
ylabel('Temperature [^{\circ} C]','FontSize',12)
title('Subsurface temperature CESM','FontSize',12)
legend({'Temperature','Polynya period'},'Location','southwest','FontSize',12)
