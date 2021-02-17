function [s,b]=qam(M,labeling)
% This function creates a square QAM constellation with M points.
% The  variable 'labeling' can be the 'BRGC' or the 'NBC'
% 
% M: Number of constellation points
% s: Constellation (Mx2)
% b: Binary labeling (Mxlog2(M))
%
% E. Agrell, M. Secondini, A. Alvarado and T. Yoshida
% Feb. 2021

%% Constellation
m=log2(M);                              % Bits per symbol
pam=[-(sqrt(M)-1):2:sqrt(M)-1].';     	% PAM constellation
s=[reshape(repmat(pam.',sqrt(M),1),1,M);reshape(repmat(pam,sqrt(M),1),1,M)].';
Es=mean(sum(abs(s).^2,2));
s=s/sqrt(Es); % Normalize to unit energy
Es=mean(sum(abs(s).^2,2));
%% Labeling
bpam=Get_Labeling(m/2,labeling);   % Binary labeling: 'BRGC' or 'NBC'
L=bin2dec(num2str(bpam));        % Convert to decimal
b1=reshape(repmat(L.',sqrt(M),1),1,M).'; % Labeling for I
b2=reshape(repmat(L,sqrt(M),1),1,M).';   % Labeling for Q
b=[dec2bin(b1,m/2)-48,dec2bin(b2,m/2)-48];% Final labeling
%figure;axis square;grid on;hold on; % To plot constellation and labeling
%for i=1:M,plot(s(i,1),s(i,2),'bx');text(s(i,1),s(i,2),num2str(b(i,:)));end

return
