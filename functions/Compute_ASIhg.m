function ASI = Compute_ASIhg(s,Pr_s,b,i,y,q)
%function ASI = Compute_ASIhg(s,Ps,b,i,y,q)
%   Compute the asymmeric information (ASI) (10) with an arbitrary
%   constellation, uniform or probabilistically-shaped symbols, 
%   soft bit-wise decoding, and an auxiliary AWGN channel. 
%   Slow, tutorial version (for-loop implementation). Histogram-based.
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
%   - ASI    :    asymmetric information, eq. (11) of the paper
%
% E. Agrell, M. Secondini, A. Alvarado and T. Yoshida
% Feb. 2021

%% Internal parameters
Bin = 32; % B in (10). can be 64 etc.
minDL = 1/32; % parameter for optimum \Delta L search (minimum value). can be 1/32, 1/16 etc.
DLmagn = 2; % parameter for optimum \Delta L search (magnification). can be 2, sqrt(2), etc.
DLcase = 5; % parameter for optimum \Delta L search (number of cases). can be 5, 10 etc.

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

%% Compute the L-value (9)
L=zeros(m,N);
aL=zeros(m,N);
for n=1:N
    for k=1:m
        num=0;
        den=0;
        for j=1:M
            if b(j,k)==0
                num=Pr_s(j)*q(y(:,n),s(:,j))+num;
            elseif b(j,k)==1
                den=Pr_s(j)*q(y(:,n),s(:,j))+den;
            end
        end
        L(k,n)=log(num./den);
        aL(k,n)=L(k,n).*((-1).^b(i(n),k));
    end
end

%% Compute the ASI (11)
aLseq = reshape(aL,size(aL,1)*size(aL,2),1);
ASItmp = zeros(DLcase,1);
for DLid = 1:1:DLcase
    scale = 1 / (minDL * (DLmagn^(DLid-1)));
    aLsq = scale_quant_L(aLseq,scale,Bin);
    ASItmp(DLid,1) = ASI_histogram_recipe(aLsq,Bin);
end
%% Choose the optimum one
id_max = find(ASItmp == max(ASItmp));
id_max1 = id_max(1);
ASI = ASItmp(id_max1,1);
return