clear all, clc

% Daan Boot, IMAU, Utrecht University
% Script to analyse data output (non-varying parameters)
% Used to make Figures 5a-f, 7a-c, 8a-c, and 9a-c

%% Determine run type
run='0';                                                % Run number
N='100';                                                % Length of simulation
n=str2num(N);                                           % Length of simulation
pm="none";                                              % Parameter
type='AF';                                              % Type run
nn=1;

%% Load in data
load(strcat(type,pm,N,run,num2str(nn)));    % Load in data
TT=75;                                      % First year
TT2=TT+25;                                  % Last year
T=365*TT+1;
T2=365*(TT2);

T_1=T_1(T:T2);                              % Select correct time period (years 76-100)
S_1=S_1(T:T2);
d_m=d_m(T:T2);
T_2=T_2(T:T2);
S_2=S_2(T:T2);

t=linspace(76,100,length(T_1));             % Time vector

%% Determine titles for plots (also represent different cases)
if isequal(type,'DF') && isequal(nn,1)
    Title='(MKL)';
elseif isequal(type,'DF') && isequal(nn,2)
    Title='(MKH)';
elseif isequal(type,'AF')
    Title='(PFB)';
elseif isequal(type,'AT')
    Title='(PFH)';
elseif isequal(type,'AS')
    Title='(PFS)';
end

%% Determine polynya periods
gT_1=gradient(T_1);
X=find(gT_1>0.15);      % Value of 0.15 is empirically determined
j=diff(X);
x=find(j==1);          
X(x)=[];

%% Plotting 5c,f and 7c, 8c, 9c
%% Create grid + density contour lines
alfa=5.82e-05;                      % Thermal expansion
beta=8e-04;                         % Haline contraction
rho0=1000;                          % Reference density

S=linspace(34.35,34.8,50);          % Salinity grid
S=repmat(S,length(S),1);

T=linspace(-2,1.3,50);              % Temperature grid
T=repmat(T,length(S),1);    

R=rho0+rho0*(-alfa*T'+beta*S);      % Density contour lines    

FS=9;

z = zeros(size(S_1));
col = linspace(TT,TT2,length(S_1)); 

gcf=figure(1)
hold on
[C,h]=contour(S(1,:),T(1,:),R,[1027.5 1027.6 1027.7 1027.8 1027.9],'Color','Black','LineWidth',0.5,'ShowText','on');
surface([S_1;S_1],[T_1;T_1],[z;z],[col;col],...
        'facecol','no',...
        'edgecol','interp',...
        'linew',2);
    
clabel(C,h,'LabelSpacing',125,'FontSize',FS)
caxis([TT TT2]);
ylim([-2 1.3])
xlim([34.35 34.8])
box on
grid on
xlabel('Salinity [g/kg]','FontSize',FS)
ylabel('Temperature [^{\circ} C]','FontSize',FS)
title(strcat('T-S',{' '},'plot',{' '},Title),'FontSize',FS)
c = colorbar;
c.Label.String = 'Years';
c.Label.FontSize=FS;
ax = gca;
ax.FontSize = FS; 

x0=500;
y0=300;
width=208;
height=150;
set(gcf,'position',[x0,y0,width,height])

%% Plot figures 5a,d, 7a, 8a, 9a
%% Density
r_1=rho0*(-alfa*T_1+beta*S_1)+rho0;
r_2=rho0*(-alfa*T_2+beta*S_2)+rho0;

%% For cases MKL and MKH below is necessary to get the correct plots

for i=1:length(r_1)
    if r_1(i)>r_2(i)
        r_2(i)=r_1(i);
        S_1(i)=S_2(i);
        T_1(i)=T_2(i);
    end
end

%% Plot
gcf=figure(2)
hold on
plot(t,r_1,'LineWidth',1.5)
plot(t,r_2,'LineWidth',1.5)

if isequal(type,'DF') && isequal(nn,1)
    Title='(MKL)';
elseif isequal(type,'DF') && isequal(nn,2)
p1=patch([t(1:end),fliplr(t(1:end))],[(1027.355*ones(1,length(t(1:end)))),fliplr((1027.9*ones(1,length(t(1:end)))))],'black','EdgeColor','k','FaceAlpha',0.1)

elseif isequal(type,'AF')
p1=patch([t(1:X(9)),fliplr(t(1:X(9)))],[(1027.355*ones(1,length(t(1:X(9))))),fliplr((1027.98*ones(1,length(t(1:X(9))))))],'black','EdgeColor','k','FaceAlpha',0.1)
p2=patch([t(X(10):X(13)),fliplr(t(X(10):X(13)))],[(1027.355*ones(1,length(t(X(10):X(13))))),fliplr((1027.98*ones(1,length(t(X(10):X(13))))))],'black','EdgeColor','k','FaceAlpha',0.1)
p3=patch([t(X(14):end),fliplr(t(X(14):end))],[(1027.355*ones(1,length(t(X(14):end)))),fliplr((1027.98*ones(1,length(t(X(14):end)))))],'black','EdgeColor','k','FaceAlpha',0.1)

elseif isequal(type,'AT')
p1=patch([t(X(1):X(end)),fliplr(t(X(1):X(end)))],[(1027.355*ones(1,length(t(X(1):X(end))))),fliplr((1027.98*ones(1,length(t(X(1):X(end))))))],'black','EdgeColor','k','FaceAlpha',0.1)

elseif isequal(type,'AS')
p1=patch([t(1:2859),fliplr(t(1:2859))],[(1027.355*ones(1,length(t(1:2859)))),fliplr((1027.98*ones(1,length(t(1:2859)))))],'black','EdgeColor','k','FaceAlpha',0.1)
p3=patch([t(X(11):end),fliplr(t(X(11):end))],[(1027.55*ones(1,length(t(X(11):end)))),fliplr((1027.98*ones(1,length(t(X(11):end)))))],'black','EdgeColor','k','FaceAlpha',0.1)    
end

grid on
box on
xlim([76 100])
ylim([1027.55 1027.8])
xlabel('Time [years]','FontSize',FS)
ylabel('Density [kgm^{-3}]','FontSize',FS)
title(strcat('\rho_1 and \rho_2',{' '},Title),'FontSize',FS)
legend({'\rho_1','\rho_2'},'Location','south','Orientation','horizontal','FontSize',FS)
ax = gca;
ax.FontSize = FS; 

xx=208;
yy=150;

x0=500;
y0=300;
width=xx;
height=yy;
set(gcf,'position',[x0,y0,width,height])

%% Plot figures 5b,e, 7b, 8b, 9b
gcf=figure(3)
hold on
if isequal(type,'DF') && isequal(nn,1)
    Title='(MKL)';
elseif isequal(type,'DF') && isequal(nn,2)
p1=patch([t(1:end),fliplr(t(1:end))],[(-1*ones(1,length(t(1:end)))),fliplr((0.95*ones(1,length(t(1:end)))))],'black','EdgeColor','k','FaceAlpha',0.1)

elseif isequal(type,'AF')
p1=patch([t(1:X(9)),fliplr(t(1:X(9)))],[(-1*ones(1,length(t(1:X(9))))),fliplr((0.95*ones(1,length(t(1:X(9))))))],'black','EdgeColor','k','FaceAlpha',0.1)
p2=patch([t(X(10):X(13)),fliplr(t(X(10):X(13)))],[(-1*ones(1,length(t(X(10):X(13))))),fliplr((0.95*ones(1,length(t(X(10):X(13))))))],'black','EdgeColor','k','FaceAlpha',0.1)
p3=patch([t(X(14):end),fliplr(t(X(14):end))],[(-1*ones(1,length(t(X(14):end)))),fliplr((0.95*ones(1,length(t(X(14):end)))))],'black','EdgeColor','k','FaceAlpha',0.1)

elseif isequal(type,'AT')
p1=patch([t(X(1):X(end)),fliplr(t(X(1):X(end)))],[(-1*ones(1,length(t(X(1):X(end))))),fliplr((0.95*ones(1,length(t(X(1):X(end))))))],'black','EdgeColor','k','FaceAlpha',0.1)

elseif isequal(type,'AS')
p1=patch([t(1:2859),fliplr(t(1:2859))],[(-1*ones(1,length(t(1:2859)))),fliplr((0.95*ones(1,length(t(1:2859)))))],'black','EdgeColor','k','FaceAlpha',0.1)
p3=patch([t(X(11):end),fliplr(t(X(11):end))],[(-1*ones(1,length(t(X(11):end)))),fliplr((0.95*ones(1,length(t(X(11):end)))))],'black','EdgeColor','k','FaceAlpha',0.1)    
end

plot(t,d_m,'LineWidth',1.5)
grid on
box on
xlim([76 100])
ylim([0 0.85])
xlabel('Time [years]','FontSize',FS)
ylabel('Sea ice thickness [m]','FontSize',FS)
title(strcat('\delta',{' '},Title),'FontSize',FS)
ax = gca;
ax.FontSize = FS; 

x0=500;
y0=300;
width=208;
height=150;
set(gcf,'position',[x0,y0,width,height])