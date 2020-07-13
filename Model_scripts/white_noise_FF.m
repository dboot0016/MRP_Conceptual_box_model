function [FF] = white_noise_FF(type,n,nn)
% Daan Boot, IMAU, Utrecht University
% This function is called on in Model.m in case of a white noise ('wn') run

F_NP=load('fnp.txt');                                       % Non-polynya freshwater flux
F_NP=repmat(F_NP,1,n);                                      % Repeat n-years
F_P=load('fp.txt');                                         % Polynya freshwater flux
F_P=repmat(F_P,1,n);                                        % Repeat n-years

if isequal(type,'DF')
    SNR_dB=4.0702;
    F_P=F_NP;
else SNR_dB=4.0702;                                         % Signal to nois ratio from CESM
end

x=F_NP;                             
W_NP=add_awgn_noise(x,SNR_dB);                              % Add white noise to non-polynya freshwater flux                        

x=F_P;
W_P=add_awgn_noise(x,SNR_dB);                               % Add white noise to polynya freshwater flux

w_NP=W_NP-(F_NP);                                           % Non-polynya white noise (with mean 0)
w_P=W_P-(F_P);                                              % Polynya white noise (with mean 0)

FF=[W_NP; W_P];

save(strcat('w_NP',type,num2str(n),num2str(nn)),'w_NP');    % Save non-polynya freshwater flux
save(strcat('w_P',type,num2str(n),num2str(nn)),'w_P');      % Save polynya noise freshwater flux

end

