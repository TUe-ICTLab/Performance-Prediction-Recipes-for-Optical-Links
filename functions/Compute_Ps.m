function Ps = Compute_Ps(s,i,y)
%function Ps = Compute_Ps(s,i,y)
%   Compute the symbol error probability Ps (3) with an arbitrary
%   constellation, uniform symbols, and a Euclidean-distance metric.
%   Tutorial version (for-loop implementation).
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
%   D,M, and N are obtained from the size of the input arguments.
% 
%   Input:
%   - s(D,M):     symbols of the M-ary modulation alphabet (M D-dimensional symbols)
%   - i(N):       sequence of transmitted symbols (indices)
%   - y(D,N):     corresponding sequence of received symbols (N D-dimensional symbols)
%   
%   Output:
%   - Ps:         symbol error rate, eq. (3) of the paper
%
% E. Agrell, M. Secondini, A. Alvarado and T. Yoshida
% Feb. 2021

%% Find and check the size of the input arguments
[D,M]=size(s);     %D=number of dimensions;M=size of the alphabet
N=length(i);    %Length of input and output sequences
[Dy,Ny]=size(y);
if (Dy~=D)
    error('Error! Output symbols and modulation symbols have different size D');
end
if (Ny~=N)
    error('Error! The input and output sequences have different length N');
end


%% Estimated index (2)
ihat=zeros(1,N);
for n=1:N
    [~,ihat(n)]=min(sum(abs(y(:,n)-s).^2,1));
end

%% SER (3)
sum_tmp=0;
for n=1:N
    sum_tmp=(ihat(n)~=i(n))+sum_tmp;
end
Ps=1/N*sum_tmp;

end

