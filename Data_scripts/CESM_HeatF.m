clear all,clc

% Daan Boot, IMAU, Utrecht University
% Script used to analyse heat flux at the atmospheric interface of CESM ru
% End product(s): Dataset for Qoa;
%                 Dataset for Qia.
% CESM files available upon request

%% Substracting variables CESM simulation
SHF=ncread('SHF_Polynya.nc','SHF');                     % Surface heat flux
If=ncread('Ice_fraction_thickness.nc','aice');          % Ice fraction    
Ih=ncread('Ice_fraction_thickness.nc','hi');            % Ice thickness
tt=ncread('Ice_fraction_thickness.nc','time');
If=If(962:2173);                                        % Select time period
Ih=Ih(962:2173);                                        % Select time period
SHF=SHF(602:1813);                                      % Select time period

%% Sort by month
Ifm=reshape(If,12,length(If)/12);                       % Sort fraction by month
Ihm=reshape(Ih,12,length(If)/12);                       % Sort thickness by month
SHFm=reshape(SHF,12,length(SHF)/12);                    % Sort heat flux by month

%% Determine condition separation Qoa and Qia (ice fraction < 0.5 --> Qoa)
X=(Ifm<50);                                             % If the fraction is lower than 0.5, the heat flux is considered to be between the ocean and the atmosphere
Y=(Ifm>=50);                                            % If the fraction is larger than 0.5, the heat flux is considered to be between the ice and the atmosphere

QOA=SHFm.*X;                                            % Use X, and Y to determine sets
QIA=SHFm.*Y;                                            % Use X, and Y to determine sets

%% Determine data sets
for i=1:12
qo=QOA(i,:);
x=find(qo==0);
qo(x)=[];
Qoa(i)=mean(qo);

qi=QIA(i,:);
y=find(qi==0);
qi(y)=[];
Qia(i)=mean(qi);
end

y=find(isnan(Qia));
Qia(y)=Qoa(y);

%% Plot data sets
Qw=[130.2 61.4 -23.6 -90.8 -143.6 -174.4 -221.8 -209.2 -84.7 -39.3 68.6 135.8];
Qi=[125.6 61.4 -23.6 -90.8 -138.8 -113.3 -94.9 -93.8 -84.7 -45.4 6.1 83.8];

xx=470;
yy=340;
FS=9;

month_vector={'Jan';'Feb';'Mar';'Apr';'May';'Jun';'Jul';'Aug';'Sep';'Oct';'Nov';'Dec'};

figure(1)
hold on
plot(1:12,Qoa,'LineWidth',1.5)
plot(1:12,Qia,'LineWidth',1.5)
%plot(1:12,Qw,'LineWidth',1.5)
%plot(1:12,Qi,'LineWidth',1.5)
grid on
box on
xlim([1 12])
xlabel('Months','FontSize',FS)
ylabel('Heat flux [Wm^{-2}]','FontSize',FS)
title('Monthly heat fluxes','FontSize',FS)
legend({'Q_{oa}','Q_{ia}'},'Orientation','Horizontal','Location','north','FontSize',FS);
ax = gca;
ax.FontSize = FS; 
set(gca,'xtick',[1:12],'xticklabel',month_vector)

x0=500;
y0=300;
width=xx;
height=yy;
set(gcf,'position',[x0,y0,width,height])
print(gcf,'appendix2.pdf','-dpdf','-r1200');

%% Interpolate monthly to daily
d=[31 28 31 30 31 30 31 31 30 31 30 31];                % Days per month

i=1;                                                    % Initial value months
for j=1:d(i)                                            % Linear interpolate two monthly values
    X1(j)=Qoa(i)+(Qoa(i+1)-Qoa(i))/d(i)*j;              % Steps are determined w.r.t. number of days in month
    Y1(j)=Qia(i)+(Qia(i+1)-Qia(i))/d(i)*j;              % For each month for both ocean-atmoshpere and ice-atmosphere sets    
end
i=i+1;
for j=1:d(i)                                            
    X2(j)=Qoa(i)+(Qoa(i+1)-Qoa(i))/d(i)*j;              
    Y2(j)=Qia(i)+(Qia(i+1)-Qia(i))/d(i)*j;              
end
i=i+1;
for j=1:d(i)                                            
    X3(j)=Qoa(i)+(Qoa(i+1)-Qoa(i))/d(i)*j;              
    Y3(j)=Qia(i)+(Qia(i+1)-Qia(i))/d(i)*j;              
end
i=i+1;
for j=1:d(i)                                            
    X4(j)=Qoa(i)+(Qoa(i+1)-Qoa(i))/d(i)*j;              
    Y4(j)=Qia(i)+(Qia(i+1)-Qia(i))/d(i)*j;              
end
i=i+1;
for j=1:d(i)                                            
    X5(j)=Qoa(i)+(Qoa(i+1)-Qoa(i))/d(i)*j;              
    Y5(j)=Qia(i)+(Qia(i+1)-Qia(i))/d(i)*j;              
end
i=i+1;
for j=1:d(i)                                            
    X6(j)=Qoa(i)+(Qoa(i+1)-Qoa(i))/d(i)*j;              
    Y6(j)=Qia(i)+(Qia(i+1)-Qia(i))/d(i)*j;              
end
i=i+1;
for j=1:d(i)                                            
    X7(j)=Qoa(i)+(Qoa(i+1)-Qoa(i))/d(i)*j;              
    Y7(j)=Qia(i)+(Qia(i+1)-Qia(i))/d(i)*j;              
end
i=i+1;
for j=1:d(i)                                            
    X8(j)=Qoa(i)+(Qoa(i+1)-Qoa(i))/d(i)*j;              
    Y8(j)=Qia(i)+(Qia(i+1)-Qia(i))/d(i)*j;              
end
i=i+1;
for j=1:d(i)                                            
    X9(j)=Qoa(i)+(Qoa(i+1)-Qoa(i))/d(i)*j;              
    Y9(j)=Qia(i)+(Qia(i+1)-Qia(i))/d(i)*j;              
end
i=i+1;
for j=1:d(i)                                            
    X10(j)=Qoa(i)+(Qoa(i+1)-Qoa(i))/d(i)*j;              
    Y10(j)=Qia(i)+(Qia(i+1)-Qia(i))/d(i)*j;              
end
i=i+1;
for j=1:d(i)                                            
    X11(j)=Qoa(i)+(Qoa(i+1)-Qoa(i))/d(i)*j;              
    Y11(j)=Qia(i)+(Qia(i+1)-Qia(i))/d(i)*j;              
end
i=i+1;
for j=1:d(i)
    X12(j)=Qoa(i)+(Qoa(1)-Qoa(i))/d(i)*j;
    Y12(j)=Qia(i)+(Qia(1)-Qia(i))/d(i)*j;
end

%%
X=[X1 X2 X3 X4 X5 X6 X7 X8 X9 X10 X11 X12];                 % A daily dataset ocean-atmosphere
Y=[Y1 Y2 Y3 Y4 Y5 Y6 Y7 Y8 Y9 Y10 Y11 Y12];                 % A daily dataset ice-atmosphere

qOA=circshift(X,15);                                        % Monthly values at 15th of each month
qIA=circshift(Y,15);                                        % To start at 1/1, shift 15 days

figure(2)
hold on
plot(1:365,qOA,'LineWidth',1.5)
plot(1:365,qIA,'LineWidth',1.5)
grid on
box on
xlim([1 365])
xlabel('Months','FontSize',12)
ylabel('Heat flux [Wm^{-2}]','FontSize',12)
title('Daily heat flux','FontSize',12)
legend({'Q_{oa}','Q_{ia}'},'Orientation','Horizontal','Location','north','FontSize',10);

%%
save('Qoa.txt','qOA','-ascii'); % The dataset for ocean-atmospheric heat fluxes consisting of daily data
save('Qia.txt','qIA','-ascii'); % The dataset for ice-atmospheric heat fluxes consisting of daily data