function Example_PhaseNoise_PS
% This code generates part of the data used for Fig. 2 of the paper
% "Performance Prediction Recipes for Optical Links", submitted to Photonics 
% Technology Letters, 2021, by Agrell, Secondini, Alvarado and Yoshida.
%
% Slow, tutorial version (for-loop implementation)
%
% Here the symbols are probabilistically shaped. For the same example, but 
% using uniform signaling, see Example_PhaseNoise_uniform.m
%
% E. Agrell, M. Secondini, A. Alvarado and T. Yoshida
% Feb. 2021

close all
addpath(genpath('functions/'))
N=1e4;                      % Number of samples
M=256;                      % Cardinality of constellation
D=2;                        % Dimensionality
m=log2(M);                  % bits per symbol
[s,b]=qam(M,'BRGC');        % Get constellation and labeling for QAM
rho=4.5923;                 % Parameter of M-B dist. Set 0 in uniform case. 4.5923 for 256QAM: H(s)=6.3 bit/2d-symbol
Pr_s=Get_MBPMF(s,rho);      % Probability mass function of s
%% Channel
sigmap2=[0.01];             % Phase noise variance
SNRdB=[5:2.5:30];           % SNR definition (only AWGN)
SNRlin=10.^(SNRdB/10);      % Linear SNR
Es=Pr_s*(sum(abs(s).^2,2)); % Average symbol energy
s=s/sqrt(Es);               % Normalization of constellation
Es=Pr_s*(sum(abs(s).^2,2)); % Average symbol energy (expected to be 1)
NoiseVar=Es./SNRlin;        % Noise variance: SNR=Es/NoiseVariance
for pp=1:length(SNRdB)
    %% TX symbols
    fprintf('SNR = %f dB \n',SNRdB(pp));
    ivec=discreteinvrnd(Pr_s,1,N); % Indices for Tx
    x=s(ivec,:).';
    %% Phase noise channel
    y=channel_phase_noise(x,NoiseVar(pp),sigmap2);
    %% Channel estimate and decoding metric
    % AWGN decoding metric
    Metric.sigma2(pp)=channel_estimate_awgn(y,x);
    qhandle_awgn=@(y,x) q_awgn(y,x,Metric.sigma2(pp));   
    % BLT decoding metric
    % The following line assumes exact knowledge of the channel parameters
    sigma2_hat(pp)=NoiseVar(pp);
    sigmap2_hat(pp)=sigmap2;
    qhandle_blt=@(y,x) q_BLT(y,x,sigma2_hat(pp),sigmap2_hat(pp)); 
    %% ASI (11) with AWGN decoding metric
    Metric.ASI(pp) = Compute_ASIhg(s.',Pr_s,b,ivec,y,qhandle_awgn);
    %% AIR_b^ps (12) with AWGN decoding metric
    Metric.AIRps(pp) = max(entropy(Pr_s) - (1-Metric.ASI(pp))*m,0);
    %% P_b^ps (13) with AWGN decoding metric
    Metric.Pbps(pp) = Compute_Pbps(s.',Pr_s,b,ivec,y,qhandle_awgn);
    %% ASI (11) with BLT decoding metric
    Metric.ASI_BLT(pp) = Compute_ASIhg(s.',Pr_s,b,ivec,y,qhandle_blt);
    %% AIR_b^ps (12) with BLT decoding metric
    Metric.AIRps_BLT(pp) = max(entropy(Pr_s) - (1-Metric.ASI_BLT(pp))*m,0);
    %% BER^ps (13) with BLT decoding metric
    Metric.Pbps_BLT(pp) = Compute_Pbps(s.',Pr_s,b,ivec,y,qhandle_blt);
end

%% Figure
figure(1)
subplot(2,1,1);
semilogy(SNRdB,Metric.Pbps,SNRdB,(1-Metric.ASI),SNRdB,Metric.Pbps_BLT,SNRdB,(1-Metric.ASI_BLT),'Linewidth',2);grid on;hold on;
axis([SNRdB(1),SNRdB(end),1e-3,1]);
hh=legend('P_b^{ps}','1-ASI','P_b^{ps} BLT','1-ASI BLT');set(hh,'Location','southwest');
title('Error probabilities');xlabel('SNR [dB]'); 

subplot(2,1,2);
plot(SNRdB,Metric.AIRps,SNRdB,Metric.AIRps_BLT,'Linewidth',2);grid on;hold on;
axis([SNRdB(1),SNRdB(end),0.8,6.2]);
hh=legend('AIR_b^{PS}','AIR_b^{PS} BLT');set(hh,'Location','northwest');
title('Achievale information rates');xlabel('SNR [dB]');

return