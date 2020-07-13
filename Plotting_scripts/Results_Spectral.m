clear all,clc

% Daan Boot, IMAU, Utrecht University
% Script to use spectral analysis on results
% Used to make Figure 10a, b, c

%% Determine run
n=100;              % Length of simulation
type='AF';          % Type of simulation (AF: Subsurface forcing; DF: no subsurface forcing)
pm='none';          % Basic parameter

%% Determine variables for analysis
X=n*365-(25*365);       % Length of analysed period (total simulation - first 25 years)
Y=100;                  % Number of different runs 
fs=1/86400;             % Sample frequency (1/day)
f =(0:X-1)*(fs/X);      % Frequency range in seconds
f=(1./(f*86400*365));   % Frequency range transformed to years
ttt=1./f;

%% Empty matrices
power_T1_wn=zeros(Y,X);
power_S1_wn=zeros(Y,X);
power_dm_wn=zeros(Y,X);
power_T2_wn=zeros(Y,X);
power_S2_wn=zeros(Y,X);

power_T1_rn=zeros(Y,X);
power_S1_rn=zeros(Y,X);
power_dm_rn=zeros(Y,X);
power_T2_rn=zeros(Y,X);
power_S2_rn=zeros(Y,X);

%% Determine power per frequency per run
for nn=1:Y;

    run='wn';                   % White noise
    load(strcat(type,pm,num2str(n),run,num2str(nn)));   % Load in data
    
    T_1=T_1(25*365+1:end);      % Skip first 25 years
    S_1=S_1(25*365+1:end);      % Skip first 25 years
    d_m=d_m(25*365+1:end);      % Skip first 25 years
    T_2=T_2(25*365+1:end);      % Skip first 25 years
    S_2=S_2(25*365+1:end);      % Skip first 25 years
    
    y_T1_wn=fft(T_1);           % Fast Fourier Transform
    y_S1_wn=fft(S_1);           % Fast Fourier Transform
    y_dm_wn=fft(d_m);           % Fast Fourier Transform
    y_T2_wn=fft(T_2);           % Fast Fourier Transform
    y_S2_wn=fft(S_2);           % Fast Fourier Transform

    power_T1_wn(nn,:)=abs(y_T1_wn).^2/X;    % Determine power
    power_S1_wn(nn,:)=abs(y_S1_wn).^2/X;    % Determine power
    power_dm_wn(nn,:)=abs(y_dm_wn).^2/X;    % Determine power
    power_T2_wn(nn,:)=abs(y_T2_wn).^2/X;    % Determine power
    power_S2_wn(nn,:)=abs(y_S2_wn).^2/X;    % Determine power
    
    T1test(nn,:)=T_1;
    
end

%% Determine mean values
% White noise
mean_T1_wn=mean(power_T1_wn,1);             
mean_S1_wn=mean(power_S1_wn,1);
mean_dm_wn=mean(power_dm_wn,1);
mean_T2_wn=mean(power_T2_wn,1);
mean_S2_wn=mean(power_S2_wn,1);

%% Determine percentiles white noise
sort_T1_wn=sort(power_T1_wn(:,2:end),1);        % Sort spectrum
T1_wn_10=sort_T1_wn(floor(Y/10)+1,:);           % 10th percentile
T1_wn_50=sort_T1_wn(floor(Y/2)+1,:);            % Median
T1_wn_90=sort_T1_wn(floor(Y*9/10)+1,:);         % 90th percentile

sort_S1_wn=sort(power_S1_wn(:,2:end),1);        % Sort spectrum
S1_wn_10=sort_S1_wn(floor(Y/10)+1,:);           % 10th percentile
S1_wn_50=sort_S1_wn(floor(Y/2)+1,:);            % Median
S1_wn_90=sort_S1_wn(floor(Y*9/10)+1,:);         % 90th percentile

sort_dm_wn=sort(power_dm_wn(:,2:end),1);        % Sort spectrum
dm_wn_10=sort_dm_wn(floor(Y/10)+1,:);           % 10th percentile
dm_wn_50=sort_dm_wn(floor(Y/2)+1,:);            % Median
dm_wn_90=sort_dm_wn(floor(Y*9/10)+1,:);         % 90th percentile

%% Plot percentiles white noise
gcf=figure(1)
hold on
p1=patch([flipud(ttt(2:end)),fliplr(flipud(ttt(2:end)))],[T1_wn_10,fliplr(T1_wn_90)],'red','EdgeColor','r','FaceAlpha',0.3);
plot(flipud(ttt(2:end)),mean_T1_wn(2:end),'LineWidth',1.5)
plot(flipud(ttt(2:end)),T1_wn_50,'k','LineWidth',1.5)
plot(flipud(ttt(2:end)),power_T1_wn(29,2:end),'Color',[0 0.5 0],'LineWidth',1.5)
plot(flipud(ttt(2:end)),T1_wn_10,'r','LineWidth',1.5)
plot(flipud(ttt(2:end)),T1_wn_90,'r','LineWidth',1.5)
set(gca,'Xscale','log')
set(gca,'Yscale','log')
xlim([ttt(2) 0.9])
ylim([10^-3 10^4])

xlabel('Frequency [years^{-1}]','FontSize',8)
ylabel('Power','FontSize',8)

x0=375;
y0=250;
width=198;
height=150;
set(gcf,'position',[x0,y0,width,height])
ax = gca;
ax.FontSize = 8; 

ax1 = gca; 
ax1_pos = ax1.Position;
ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','left',...
    'Color','none','XTick',[2,3,5,10,20,30,50],'YTick',[]);

line(fliplr(f(2:end)),mean_T1_wn(2:end),'Parent',ax2)

xlim([1/0.9 1/ttt(2)])
ylim([10^-3 10^4])
set(gca,'Yscale','log')
set ( gca, 'xdir', 'reverse' )
set(gca,'Xscale','log')
hold off
grid on
box on
xlabel('Period [years]','FontSize',8)
ax = gca;
ax.FontSize = 8; 

title('Spectral analysis \delta (extra: n=100)','FontSize',8)

x0=375;
y0=250;
width=198;
height=150;
set(gcf,'position',[x0,y0,width,height])

%%
gcf=figure(2)
hold on
p1=patch([flipud(ttt(2:end)),fliplr(flipud(ttt(2:end)))],[S1_wn_10,fliplr(S1_wn_90)],'red','EdgeColor','r','FaceAlpha',0.3);
plot(flipud(ttt(2:end)),mean_S1_wn(2:end),'LineWidth',1.5)
plot(flipud(ttt(2:end)),S1_wn_50,'k','LineWidth',1.5)
plot(flipud(ttt(2:end)),power_S1_wn(29,2:end),'Color',[0 0.5 0],'LineWidth',1.5)
plot(flipud(ttt(2:end)),S1_wn_10,'r','LineWidth',1.5)
plot(flipud(ttt(2:end)),S1_wn_90,'r','LineWidth',1.5)
set(gca,'Xscale','log')
set(gca,'Yscale','log')
xlim([ttt(2) 0.9])
ylim([10^-4 10^2])

xlabel('Frequency [years^{-1}]','FontSize',8)
ylabel('Power','FontSize',8)

x0=375;
y0=250;
width=198;
height=150;
set(gcf,'position',[x0,y0,width,height])
ax = gca;
ax.FontSize = 8; 

ax1 = gca; 
ax1_pos = ax1.Position;
ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','left',...
    'Color','none','XTick',[2,3,5,10,20,30,50],'YTick',[]);

line(fliplr(f(2:end)),mean_S1_wn(2:end),'Parent',ax2)

xlim([1/0.9 1/ttt(2)])
ylim([10^-4 10^2])
set(gca,'Yscale','log')
set ( gca, 'xdir', 'reverse' )
set(gca,'Xscale','log')
hold off
grid on
box on
xlabel('Period [years]','FontSize',8)
ax = gca;
ax.FontSize = 8; 

title('Spectral analysis \delta (extra: n=100)','FontSize',8)


x0=375;
y0=250;
width=198;
height=150;
set(gcf,'position',[x0,y0,width,height])

%%
gcf=figure(3)
hold on
p1=patch([flipud(ttt(2:end)),fliplr(flipud(ttt(2:end)))],[dm_wn_10,fliplr(dm_wn_90)],'red','EdgeColor','r','FaceAlpha',0.3);
plot(flipud(ttt(2:end)),mean_dm_wn(2:end),'LineWidth',1.5)
plot(flipud(ttt(2:end)),dm_wn_50,'k','LineWidth',1.5)
plot(flipud(ttt(2:end)),power_dm_wn(29,2:end),'Color',[0 0.5 0],'LineWidth',1.5)
plot(flipud(ttt(2:end)),dm_wn_10,'r','LineWidth',1.5)
plot(flipud(ttt(2:end)),dm_wn_90,'r','LineWidth',1.5)
set(gca,'Xscale','log')
set(gca,'Yscale','log')
xlim([ttt(2) 0.9])
ylim([10^-6 10^2])

xlabel('Frequency [years^{-1}]','FontSize',8)
ylabel('Power','FontSize',8)

x0=375;
y0=250;
width=198;
height=150;
set(gcf,'position',[x0,y0,width,height])
ax = gca;
ax.FontSize = 8; 

ax1 = gca; 
ax1_pos = ax1.Position;
ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','left',...
    'Color','none','XTick',[2,3,5,10,20,30,50],'YTick',[]);

line(fliplr(f(2:end)),mean_dm_wn(2:end),'Parent',ax2)

xlim([1/0.9 1/ttt(2)])
ylim([10^-6 10^2])
set(gca,'Yscale','log')
set ( gca, 'xdir', 'reverse' )
set(gca,'Xscale','log')
hold off
grid on
box on
xlabel('Period [years]','FontSize',8)
ax = gca;
ax.FontSize = 8; 

title('Spectral analysis \delta (extra: n=100)','FontSize',8)


x0=375;
y0=250;
width=198;
height=150;
set(gcf,'position',[x0,y0,width,height])
