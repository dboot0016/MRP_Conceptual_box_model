clear all,clc

% Daan Boot, IMAU, Utrecht University
% Script used to analyse freshwater flux at the atmospheric interface of CESM ru
% End product(s): Dataset for FNP;
%                 Dataset for FP.
% CESM files available upon request

%% Substract variables from CESM simulation
t=ncread('Polynya_PREC_EVAP.nc','time');                % Time vector oceanic component
F=ncread('Polynya_PREC_EVAP.nc','NET');                 % Net precipitation vector oceanic component
P=ncread('Polynya_PREC_EVAP.nc','PREC');                % Precipitation vector oceanic component
E=ncread('Polynya_PREC_EVAP.nc','EVAP');                % Evaporation vector oceanic component

ta=ncread('Polynya_PREC_EVAP_Atm.nc','time');           % Time vector atmospheric component
Fa=ncread('Polynya_PREC_EVAP_Atm.nc','NET');            % Net precipitation vector atmospheric component
Pa=ncread('Polynya_PREC_EVAP_Atm.nc','PREC');           % Precipitation vector atmospheric component
Ea=ncread('Polynya_PREC_EVAP_Atm.nc','EVAP');           % Evaporation vector atmospheric component

F=F(602:1813);                                          % Select years 150-250
P=P(602:1813);                                          % Select years 150-250
E=E(602:1813);                                          % Select years 150-250
        
Fa=Fa(602:1813);                                        % Select years 150-250
Pa=Pa(602:1813);                                        % Select years 150-250
Ea=Ea(602:1813);                                        % Select years 150-250

%% Determining freshwater input
Polynya=[337:396 661:708 973:1032];                     % Determine polynya periods
N=1:1212;                                               % Entire time vector
N(Polynya)=[];                                          % Determine non-polynya periods

%% Determine mean monthly values: oceanic component
FNP=F(N);                                               % Select net precipitation vector for non-polynya years 
FP=F(Polynya);                                          % Select net precipitation vector for polynya years

FNP1=reshape(FNP,12,length(FNP)/12);                    % Order by month
FP1=reshape(FP,12,length(FP)/12);                       % Order by month
F1=reshape(F,12,length(F)/12);                          % Order by month (total dataset)

for i=1:12                                              % Take mean per month
    FoNP(i)=mean(FNP1(i,:));
    FoP(i)=mean(FP1(i,:));
    FoT(i)=mean(F1(i,:));
end

%% Determine mean monthly values: atmospheric component
FaNP=Fa(N);                                             % Select net precipitation vector for non-polynya years
FaP=Fa(Polynya);                                        % Select net precipitation vector for polynya years

FaNP1=reshape(FaNP,12,length(FaNP)/12);                 % Order by month
FaP1=reshape(FaP,12,length(FaP)/12);                    % Order by month
Fa1=reshape(Fa,12,length(Fa)/12);                       % Order by month

for i=1:12                                              % Take monthly mean
    Fnp(i)=mean(FaNP1(i,:));
    Fp(i)=mean(FaP1(i,:));
    Ft(i)=mean(Fa1(i,:));
end

%% Plot freshwater input: oceanic component
figure(1)
hold on
plot(1:12,FoNP,'LineWidth',1.5)
plot(1:12,FoP,'LineWidth',1.5)
plot(1:12,FoT,'LineWidth',1.5)
grid on
box on
xlim([1 12])
xlabel('Months','FontSize',12)
ylabel('Freshwater input [mm/day]','FontSize',12)
title('Montlhy freshwater input','FontSize',12)
legend({'F_{NP}','F_P','F_{total}'},'Orientation','Horizontal','Location','northeast','FontSize',10);

%% Plot freshwater input: atmospheric component
xx=470;
yy=340;
FS=9;

hoi={'Jan';'Feb';'Mar';'Apr';'May';'Jun';'Jul';'Aug';'Sep';'Oct';'Nov';'Dec'};
figure(2)
hold on
plot(1:12,Fnp,'LineWidth',1.5)
plot(1:12,Fp,'LineWidth',1.5)
%plot(1:12,Ft,'LineWidth',1.5)
grid on
box on
xlim([1 12])
xlabel('Months','FontSize',FS)
ylabel('Freshwater input [mm/day]','FontSize',FS)
title('Montlhy freshwater input','FontSize',FS)
legend({'F_{NP}','F_P','F_{total}'},'Orientation','Horizontal','Location','southwest','FontSize',10);
ax = gca;
ax.FontSize = FS; 
set(gca,'xtick',[1:12],'xticklabel',hoi)

x0=500;
y0=300;
width=xx;
height=yy;
set(gcf,'position',[x0,y0,width,height])
%print(gcf,'appendix1.pdf','-dpdf','-r1200');

%% Plot freshwater input: compare new with old
FoldNP=[1.46 1.81 2.23 2.34 2.27 1.89 1.62 1.43 1.43 1.66 1.76 1.60];      % Old dataset non-polynya years (Martinson et al., 1981)
FoldP=[1.47 1.83 2.26 2.37 2.27 1.98 1.69 1.52 1.50 1.71 1.74 1.60];       % Old dataset polynya years (Martinson et al., 1981)

figure(3)
hold on
plot(1:12,Fnp,'LineWidth',1.5)
plot(1:12,Fp,'LineWidth',1.5)
plot(1:12,FoldNP,'LineWidth',1.5)
plot(1:12,FoldP,'LineWidth',1.5)
grid on
box on
xlim([1 12])
ylim([0 3])
xlabel('Months','FontSize',12)
ylabel('Freshwater input [mm/day]','FontSize',12)
title('Montlhy freshwater input','FontSize',12)
legend({'F_{NP}','F_P','F_{NP} (old)','F_P (old)'},'Orientation','Horizontal','Location','north','FontSize',10);

%% Transform into salinity flux
fnp=Fnp*1e-03/86400*35;         % Unit transformation
fp=Fp*1e-03/86400*35;           % Unit transformation

figure(4)
hold on
plot(1:12,fnp,'LineWidth',1.5)
plot(1:12,fp,'LineWidth',1.5)
grid on
box on
xlim([1 12])
xlabel('Months','FontSize',12)
ylabel('Freshwater flux [g/kg s^{-1}]','FontSize',12)
title('Montlhy freshwater flux','FontSize',12)
legend({'F_{NP}','F_P'},'Orientation','Horizontal','Location','north','FontSize',10);

%% Interpolate monthly to daily
d=[31 28 31 30 31 30 31 31 30 31 30 31];                % Days per month

i=1;                                                    % Initial value months
for j=1:d(i)                                            % Linear interpolate two monthly values
    X1(j)=fnp(i)+(fnp(i+1)-fnp(i))/d(i)*j;              % Steps are determined w.r.t. number of days in month
    Y1(j)=fp(i)+(fp(i+1)-fp(i))/d(i)*j;                 % For each month for both non-polynya and polynya years    
end
i=i+1;
for j=1:d(i)
    X2(j)=fnp(i)+(fnp(i+1)-fnp(i))/d(i)*j;
    Y2(j)=fp(i)+(fp(i+1)-fp(i))/d(i)*j;
end
i=i+1;
for j=1:d(i)
    X3(j)=fnp(i)+(fnp(i+1)-fnp(i))/d(i)*j;
    Y3(j)=fp(i)+(fp(i+1)-fp(i))/d(i)*j;
end
i=i+1;
for j=1:d(i)
    X4(j)=fnp(i)+(fnp(i+1)-fnp(i))/d(i)*j;
    Y4(j)=fp(i)+(fp(i+1)-fp(i))/d(i)*j;
end
i=i+1;
for j=1:d(i)
    X5(j)=fnp(i)+(fnp(i+1)-fnp(i))/d(i)*j;
    Y5(j)=fp(i)+(fp(i+1)-fp(i))/d(i)*j;
end
i=i+1;
for j=1:d(i)
    X6(j)=fnp(i)+(fnp(i+1)-fnp(i))/d(i)*j;
    Y6(j)=fp(i)+(fp(i+1)-fp(i))/d(i)*j;
end
i=i+1;
for j=1:d(i)
    X7(j)=fnp(i)+(fnp(i+1)-fnp(i))/d(i)*j;
    Y7(j)=fp(i)+(fp(i+1)-fp(i))/d(i)*j;
end
i=i+1;
for j=1:d(i)
    X8(j)=fnp(i)+(fnp(i+1)-fnp(i))/d(i)*j;
    Y8(j)=fp(i)+(fp(i+1)-fp(i))/d(i)*j;
end
i=i+1;
for j=1:d(i)
    X9(j)=fnp(i)+(fnp(i+1)-fnp(i))/d(i)*j;
    Y9(j)=fp(i)+(fp(i+1)-fp(i))/d(i)*j;
end
i=i+1;
for j=1:d(i)
    X10(j)=fnp(i)+(fnp(i+1)-fnp(i))/d(i)*j;
    Y10(j)=fp(i)+(fp(i+1)-fp(i))/d(i)*j;
end
i=i+1;
for j=1:d(i)
    X11(j)=fnp(i)+(fnp(i+1)-fnp(i))/d(i)*j;
    Y11(j)=fp(i)+(fp(i+1)-fp(i))/d(i)*j;
end
i=i+1;
for j=1:d(i)
    X12(j)=fnp(i)+(fnp(1)-fnp(i))/d(i)*j;
    Y12(j)=fp(i)+(fp(1)-fp(i))/d(i)*j;
end

%%
X=[X1 X2 X3 X4 X5 X6 X7 X8 X9 X10 X11 X12];                 % A daily dataset non-polynya
Y=[Y1 Y2 Y3 Y4 Y5 Y6 Y7 Y8 Y9 Y10 Y11 Y12];                 % A daily dataset polynya

FFNP=circshift(X,15);                                       % Monthly values at 15th of each month
FFP=circshift(Y,15);                                        % To start at 1/1, shift 15 days

figure(5)
hold on
plot(1:365,FFNP,'LineWidth',1.5)
plot(1:365,FFP,'LineWidth',1.5)
grid on
box on
xlim([1 365])
xlabel('Months','FontSize',12)
ylabel('Freshwater flux [g/kg s^{-1}]','FontSize',12)
title('Daily freshwater flux','FontSize',12)
legend({'F_{NP}','F_P'},'Orientation','Horizontal','Location','north','FontSize',10);

%%
save('fnp.txt','FFNP','-ascii');    % The eventual data set with daily values for non-polynya years 
save('fp.txt','FFP','-ascii');      % The eventual data set with daily values for polynya years
