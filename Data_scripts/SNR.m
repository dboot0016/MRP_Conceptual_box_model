clear all, clc

% Daan Boot, IMAU, Utrecht University
% Script to determine SNR from CESM data
% CESM data available upon request

n=150;                                      % Number of years we want the noise for
t=linspace(1,n,n*365);                      % Time vector for plotting
%% Part 1: White noise

%% Load in data
F_CESM=ncread('Polynya_PREC_EVAP_Atm.nc','NET');            % Net precipitation vector atmospheric component
F_CESM=F_CESM(602:1813);                                    % Select years 150-250

F_NP=load('fnp.txt');                                       % Non-polynya freshwater flux
F_NP=repmat(F_NP,1,n);                                      % Repeat n-years
F_P=load('fp.txt');                                         % Polynya freshwater flux
F_P=repmat(F_P,1,n);                                        % Repeat n-years

%% Determine SN-ratio
SNR_dB=mean(F_CESM.^2)/var(F_CESM);                         % Signal to noise ratio
