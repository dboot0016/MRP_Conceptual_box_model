clear all,clc

% Daan Boot, IMAU, Utrecht University
% Used to make T1-S1 plot of CESM data (Figure 6)
% CESM data available upon request

%% Data
T=ncread('TEMP_Polynya_depth.nc','TEMP');   % CESM data temperature
T=T(:,602:1813);                            % CESM data years 150-250
l1=3; l2=12;                                % Indices for first layer (depth)
T1=T(l1:l2,:);                              % Temperature layer 1

S=ncread('SALT_Polynya_depth.nc','SALT');   % CESM data salinity
S=S(:,602:1813);                            % CESM data years 150-250
l1=3; l2=12;                                % Indices for first layer (depth)
S1=S(l1:l2,:);                              % Temperature layer 1

l=ncread('Layer_depth.nc','layer');         % Size per layer [m]

t21=408; t22=708;                           % Time indices period 2

%% Average layer characteristics
T1m=sum((T1(:,:)).*l(l1:l2))/sum(l(l1:l2)); % Average layer temperature
S1m=sum((S1(:,:)).*l(l1:l2))/sum(l(l1:l2)); % Average layer salinity
T1_1=T1m(:,t21:t22);                        % Select time period
S1_1=S1m(:,t21:t22);                        % Select time period

%% Preparing vectors for axis plotting
TT=210;                                     % Starting model year
TT2=TT+25;                                  % Ending model year

x = S1_1;
y = T1_1;
z = zeros(size(S1_1));
col = linspace(TT,TT2,length(S1_1));        % This is the color, vary with x in this case.

S=linspace(min(S1_1)-0.05,max(S1_1)+0.05,50);
S=sort(S);
S=repmat(S,length(S),1);                    % Salinity vector for x-axis

T=linspace(min(T1_1)-0.2,max(T1_1)+0.2,50);
T=sort(T);
T=repmat(T,length(S),1);                    % Temperature vector for y-axis

alfa=5.82e-05;
beta=8e-04;
rho0=1000;

R=rho0+rho0*(-alfa*T'+beta*S);              % Density contour lines
%% Plotting
xx=520;
yy=380;
FS=10;

figure(3)
hold on
[C,h]=contour(S(1,:),T(1,:),R,[1027.5 1027.6 1027.7 1027.8 1027.9 1028],'Color','Black','LineWidth',0.5,'ShowText','on');
surface([x;x],[y;y],[z;z],[col;col],...
        'facecol','no',...
        'edgecol','interp',...
        'linew',2);
    
clabel(C,h,'LabelSpacing',180)
caxis([TT TT2]);
box on
grid on
xlabel('Salinity [g/kg]','FontSize',FS)
ylabel('Temperature [^{\circ} C]','FontSize',FS)
title('T-S plot (CESM [210-235])','FontSize',FS)
c = colorbar;
c.Label.String = 'CESM model years';
c.Label.FontSize=FS;
ax = gca;
ax.FontSize = FS; 

x0=500;
y0=300;
width=xx;
height=yy;
set(gcf,'position',[x0,y0,width,height])
