function var = channel_estimate_awgn(y,x)
%function var = channel_estimate_awgn(y,x)
%   Estimate the variance (per dimension) of the AWGN channel model
%   given the input and output realizations (x and y). This is (1) in 
% "Performance Prediction Recipes for Optical Links", Photonics Technology
% Letters, 2021, by Agrell, Secondini, Alvarado and Yoshida.
%
%   y and x are real DxN matrices (N D-dimensional symbols)
%   The channel model is y=x+n, where n is a DxN matrix of zero-mean
%   Gaussian variables
%   The variance is estimated by averaging the square of the noise
%   components over time and dimensions.
%
% E. Agrell, M. Secondini, A. Alvarado and T. Yoshida
% Feb. 2021

var=mean(mean(abs(y-x).^2));  %squared error, averaged over dimensions and time          

end

