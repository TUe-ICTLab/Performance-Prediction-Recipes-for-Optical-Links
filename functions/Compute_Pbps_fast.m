function Pbps = Compute_Pbps_fast(s,Pr_s,b,i,y,q)
%function Pbps = Compute_Pbps_fast(s,Pr_s,b,i,y,q)
%   Compute the pre-FEC BER (12) with an arbitrary constellation, 
%   uniform or probabilistically-shaped symbols, 
%   hard bit-wise decoding, and an auxiliary AWGN channel. 
%   Relatively fast (matrix implementation).
%
%   The input and output sequences are passed as input arguments and must 
%   be generated outside this function, experimentally or by a numerical 
%   simulation of the true channel.
%
%   The M-ary modulation alphabet is made of M real vector symbols in a
%   D-dimensional space(e.g., D=1 for real symbols, D=2 for complex symbols, 
%   D=4 for dual-polarization complex symbols, ...).
%
%   The input sequence consists of N indices, i={i1,i2,...}. Each index 
%   is mapped to one of the M symbols of the alphabet.
%
%   The output sequence consists of N real vector symbols in a
%   D-dimensional spaces.
%
%   D, M, N, and m are obtained from the size of the input arguments.
% 
%   Input:
%   - s(D,M) :    symbols of the M-ary modulation alphabet (M D-dimensional symbols)
%   - Pr_s(M):    symbol probabilities (M D-dimensional symbols)
%   - b(M,m) :    labelling of bits for symbol
%   - i(N)   :    sequence of transmitted symbols (indices)
%   - y(D,N) :    corresponding sequence of received symbols (N D-dimensional symbols)
%   - q      :    function handle for the auxiliary channel q(y,s)
%   
%   Output:
%   - Pbps   :    pre-FEC BER, eq. (13) of the paper
%
% E. Agrell, M. Secondini, A. Alvarado and T. Yoshida
% Feb. 2021

%% Find and check the size of the input arguments
[D,M]=size(s);  % D=number of dimensions;M=size of the alphabet
N=length(i);    % Length of input and output sequences
m=size(b,2);    % Number of bits per symbol
[Dy,Ny]=size(y);
if (Dy~=D)
    error('Error! Output symbols and modulation symbols have different size D');
end
if (Ny~=N)
    error('Error! The input and output sequences have different length N');
end

%% Use an AWGN model for q(y,s) if not passed as argument
if nargin<6
    warning('Warning! Auxiliary channel not defined. Using AWGN with variance estimated from given data');
    sigma2=mean(mean(abs(y-s(:,i)).^2)); 
    q=@(y,x) q_awgn(y,x,sigma2);  
end

%% Compute the L-value (9) and the pre-FEC BER (13)
sum_tmp=zeros(1,N);
aL=zeros(m,N);
for k=1:m
    num=zeros(1,N);
    den=zeros(1,N);
    for j=1:M
        if b(j,k)==0
            num=Pr_s(j)*q(y(:,:),s(:,j))+num;
        elseif b(j,k)==1
            den=Pr_s(j)*q(y(:,:),s(:,j))+den;
        end
    end
    aL(k,:)=log(num./den).*(((-1).^b(i,k)).');
    sum_tmp=sum_tmp+(aL(k,:)<=0);
end
Pbps=(1/(m*N))*sum(sum_tmp);
return