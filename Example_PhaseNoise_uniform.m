function Example_PhaseNoise_uniform
% This code generates part of the data used for Fig. 2 of the paper
% "Performance Prediction Recipes for Optical Links", Photonics Technology
% Letters, 2021, by Agrell, Secondini, Alvarado and Yoshida.
%
% Here the symbols are equally likely. For the same example, but using
% probabilistic shaping, see Example_PhaseNoise_PS.m
%
% E. Agrell, M. Secondini, A. Alvarado and T. Yoshida
% Feb. 2021

close all
addpath(genpath('functions/'))
N=1e4;                      % Number of samples
M=64;                       % Cardinality of constellation
D=2;                        % Dimensionality
m=log2(M);                  % bits per symbol
[s,b]=qam(M,'BRGC');         % Get constellation and labeling for QAM
%% Channel
sigmap2=[0.01];             % Variance of Phase noise
SNRdB=[5:2.5:30];            %  SNR definition (only AWGN)
SNRlin=10.^(SNRdB/10);      
Es=mean(sum(abs(s).^2,2));
NoiseVar=Es./SNRlin;        % AWGN Noise variance: SNR=Es/NoiseVariance
for pp=1:length(SNRdB)
    %% TX symbols
    fprintf('SNR = %f dB \n',SNRdB(pp));
    ivec=randi(M,1,N);          % Indices for Tx
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
    %% SER (3)
    Metric.Ps(pp) = Compute_Ps(s.',ivec,y);
    %% BER (4)
    Metric.Pb(pp) = Compute_Pb(s.',b,ivec,y);
    %% AIR_b with HD (last paragraph of Sec. III)
    Metric.AIRb_HD(pp)=m*(1-binary_entropy(Metric.Pb(pp)));
    %% AIR_s (5) with AWGN decoding metric
    Metric.AIRs(pp) = Compute_AIRs(s.',ivec,y,qhandle_awgn);
    %% AIR_b (6) with AWGN decoding metric
    Metric.AIRb(pp) = Compute_AIRb(s.',b,ivec,y,qhandle_awgn);
    %% AIR_s (5) with BLT decoding metric
    Metric.AIRs_BLT(pp) = Compute_AIRs(s.',ivec,y,qhandle_blt);
    %% AIR_b (6) with BLT decoding metric
    Metric.AIRb_BLT(pp) = Compute_AIRb(s.',b,ivec,y,qhandle_blt);
    %% Plot received constellation (for visualization purposes only)
    figure(1);subplot(3,ceil(length(SNRdB)/3),pp);hold on;axis square;grid on;
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    for i=1:M,pnt=find(ivec==i);plot(y(1,pnt),y(2,pnt),'o','Color',rand(3,1));end
    title(['SNR=',num2str(SNRdB(pp))]);axis([-2,2,-2,2]);
    pause(0.1)
end

figure(2);
subplot(2,1,1);title('Error probabilities');xlabel('SNR [dB]'); 
semilogy(SNRdB,Metric.Ps,SNRdB,Metric.Pb,'Linewidth',2);grid on;hold on;
axis([SNRdB(1),SNRdB(end),1e-3,1]);
hh=legend('P_s','P_b');set(hh,'Location','southwest')
subplot(2,1,2);title('Achievale information rates');xlabel('SNR [dB]');
plot(SNRdB,Metric.AIRb_HD,SNRdB,Metric.AIRs,SNRdB,Metric.AIRb,SNRdB,Metric.AIRs_BLT,SNRdB,Metric.AIRb_BLT,'Linewidth',2);
grid on;hold on;
axis([SNRdB(1),SNRdB(end),0.8,6.2]);
hh=legend('AIRb-HD','AIR_s','AIR_b','AIR_s BLT','AIR_b BLT');
set(hh,'Location','northwest')
set(gcf,'units','normalized','outerposition',[0 0 1 1]);

return
