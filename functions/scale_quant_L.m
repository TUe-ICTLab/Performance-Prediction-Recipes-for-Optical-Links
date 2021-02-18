function[Lsq] = scale_quant_L(L,scale,Bin)
% function[Lsq] = scale_quant_L(L,scale,Bin)
% normalization and quantization of L-value
%
%   Input:
%   - L     :    input L-value
%   - scale :    scaling parameter
%   - Bin   :    number of bins in L-value quantization
%   
%   Output:
%   - Lsq   :    scaled and quantized L-value
%
% E. Agrell, M. Secondini, A. Alvarado and T. Yoshida
% Feb. 2021
tmp = L*scale;
halfBin = (Bin-mod(Bin,2))/2;
if mod(Bin,2)==0 % output: +/-1,+/-2,...,+/-(Bin/2)
    Lsq = sign(tmp).*floor(abs(tmp)+1);
    Lsq = max(Lsq,-halfBin);
    Lsq = min(Lsq,halfBin);
else % output: 0,+/-1,+/-2,...,+/-((Bin-1)/2)
    Lsq = sign(tmp).*floor(abs(tmp));
    Lsq = max(Lsq,-halfBin);
    Lsq = min(Lsq,halfBin);    
end
return