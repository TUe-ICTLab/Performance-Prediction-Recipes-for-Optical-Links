function AIRb = Compute_AIRb(s,b,i,y,q)
%function AIRb = Compute_AIRb(s,b,i,y,q)
%   Compute the achievable information rate AIRb (6) with an arbitrary
%   constellation, uniform symbols, and soft bit-wise decoding.
%   Tutorial version (for-loop implementation).
%
%   The input and output sequences are passed as input arguments and must 
%   be generated outside this function, experimentally or by a numerical 
%   simulation of the true channel.
%
%   The auxiliary channel q(y,s) is passed as an input argument (function
%   handle). If absent, an AWGN auxiliary channel is used by default and
%   its variance is estimated from data.
%
%   The M-ary modulation alphabet is made of M real vector symbols in a
%   D-dimensional space(e.g., D=1 for real symbols, D=2 for complex symbols, 
%   D=4 for dual-polarization complex symbols, ...).
%
%   The binary labelling consists of M labels of m=log2M bits.
%
%   The input sequence consists of N indices, i={i1,i2,...}. Each index 
%   is mapped to one of the M symbols of the alphabet.
%
%   The output sequence consists of N real vector symbols in a
%   D-dimensional spaces.
%
%   D,M, and N are obtained from the size of the input arguments.
% 
%   The variance (per dimension) of the auxiliary AWGN channel is passed as 
%   an optional input argument or, if absent, estimated internally.
%
%   Input:
%   - s(D,M):     symbols of the M-ary modulation alphabet (M D-dimensional symbols)
%   - b(M,m):     binary labels (m bits for each of the M symbols)
%   - i(N):       sequence of transmitted symbols (indices)
%   - y(D,N):     corresponding sequence of received symbols (N D-dimensional symbols)
%   - q:          function handle for the auxiliary channel q(y,s)
%   
%   Output:
%   - AIRb:       achievable information rate (bit/D-dim symbol), eq. (6) of the paper
%
% E. Agrell, M. Secondini, A. Alvarado and T. Yoshida
% Feb. 2021

%% Find and check the size of the input arguments
[D,M]=size(s);     %D=number of dimensions;M=size of the alphabet
[Mb,m]=size(b);     %m must be log2(M) (number of bits)
N=length(i);    %Length of input and output sequences
[Dy,Ny]=size(y);
if (Dy~=D)
    error('Error! Output symbols and modulation symbols have different size D');
end
if (Ny~=N)
    error('Error! The input and output sequences have different length N');
end
if (m~=log2(M))
    error('Error! The number of bits must be m=log2(M), where M is the size of the alphabet');
end
if (Mb~=M)
    error('Error! The number of symbols and the number of binary labels should be the same (M)');
end

%% Use an AWGN model for q(y,s) if not passed as argument
if nargin<5
    warning('Warning! Auxiliary channel not defined. Using AWGN with variance estimated from given data');
    sigma2=mean(mean(abs(y-s(:,i)).^2)); 
    q=@(y,x) q_awgn(y,x,sigma2);  
end

%% Compute the AIR
sum_tmp=0;
for n=1:N
    for k=1:m
        num=0;
        den=0;
        for j=1:M
            q_y_s=q(y(:,n),s(:,j));
            num=q_y_s+num;
            if b(j,k)==b(i(n),k)
                den=q_y_s+den;
            end
        end
        sum_tmp=sum_tmp+log2(num./den);
    end
end
AIRb=m-(1/N)*sum_tmp;

return
