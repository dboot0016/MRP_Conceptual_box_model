clear all,clc

% Daan Boot, IMAU, Utrecht University
% Script to plot annual mixed layer depth and potential density
% CESM data available upon request
% Used to make Figure 2a, b

%% Load in variables
MX=ncread('Mixed_layer_depth_Polynya.nc','XMXL');           % Mixed layer depth
MX=MX(602:1813);                                            % Select years 150-250                                           

t=linspace(150,250,length(MX));                             % Time vector

%% Plotting
xx=470;
yy=340;
FS=9;
figure(1)
plot(t,-MX,'LineWidth',1.5)

grid on
box on
xlabel('Time [CESM model years]','FontSize',FS)
ylabel('Maximum mixed layer depth [m]','FontSize',FS)
title('Maximum mixed layer depth','FontSize',FS)
%legend({'F_{NP}','F_P','F_{total}'},'Orientation','Horizontal','Location','southwest','FontSize',10);
ax = gca;
ax.FontSize = FS; 
xticks([150 160 170 180 190 200 210 220 230 240 250])
yticks([-2000 -1750 -1500 -1250 -1000 -750 -500 -250 0])

yticklabels({'2000','1750','1500','1250','1000','750','500','250','0'})


x0=500;
y0=300;
width=xx;
height=yy;
set(gcf,'position',[x0,y0,width,height])
%print(gcf,'appendix3.pdf','-dpdf','-r1200');

%% Potential density
PD=ncread('POT_DENS_Polynya_depth.nc','RHO');   % Potential density           
PD=PD(1:41,602:1813);                           % Select layers 1-41, years 150-250

d=ncread('Layer_depth.nc','depth');         % Mean depth at cell
D=-d(1:41);                                 % Mean depth converted for plotting

pd=movmean(PD,60,2);                        % Take moving average of 5 years (60 months)

t=linspace(150,250,length(PD));             % Time vector

%% Plotting
xx=470;
yy=340;
FS=9;

figure(1)
contourf(t,D(1:30),pd(1:30,:))
grid on
box on
xlabel('Time [CESM model years]','FontSize',FS)
ylabel('Depth [m]','FontSize',FS)
title('Potential density','FontSize',FS)

ax = gca;
ax.FontSize = FS; 
xticks([150 160 170 180 190 200 210 220 230 240 250])
yticks([-2750 -2500 -2250 -2000 -1750 -1500 -1250 -1000 -750 -500 -250 0])

yticklabels({'2750','2500','2250','2000','1750','1500','1250','1000','750','500','250','0'})
c = colorbar;
c.Label.String = 'Potential density [kg m^{-3}]';

x0=500;
y0=300;
width=xx;
height=yy;
set(gcf,'position',[x0,y0,width,height])
%print(gcf,'appendix4.pdf','-dpdf','-r1200');