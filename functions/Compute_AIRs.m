function AIRs = Compute_AIRs(s,i,y,q)
%function AIRs = Compute_AIRs(s,i,y,q)
%   Compute the achievable information rate AIRs (5) with an arbitrary
%   constellation, uniform symbols, and soft symbol-wise decoding.
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
%   The input sequence consists of N indices, i={i1,i2,...}. Each index 
%   is mapped to one of the M symbols of the alphabet.
%
%   The output sequence consists of N real vector symbols in a
%   D-dimensional spaces.
%
%   D,M, and N are obtained from the size of the input arguments.
% 
%
%   Input:
%   - s(D,M):     symbols of the M-ary modulation alphabet (M D-dimensional symbols)
%   - i(N):       sequence of transmitted symbols (indices)
%   - y(D,N):     corresponding sequence of received symbols (N D-dimensional symbols)
%   - q:          function handle for the auxiliary channel q(y,s)
%   
%   Output:
%   - AIRs:       achievable information rate (bit/D-dim symbol), eq. (5) of the paper
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

%% Use an AWGN model for q(y,s) if not passed as argument
if nargin<4
    warning('Warning! Auxiliary channel not defined. Using AWGN with variance estimated from given data');
    sigma2=mean(mean(abs(y-s(:,i)).^2)); 
    q=@(y,x) q_awgn(y,x,sigma2);  
end

%% Compute the AIR
sum_tmp=0;
for n=1:N
    q_y_giv_x=q(y(:,n),s(:,i(n)));
    q_y=0;
    for j=1:M
        q_y=q_y+q(y(:,n),s(:,j));
    end
    sum_tmp=sum_tmp+log2(q_y_giv_x./q_y);
end
AIRs=log2(M)+sum_tmp/N;

return
