function [Xout, f] = ProtonIPFading(t, X, varargin)
% Corrects for image plate fading for three image plate types: tr (blue,
% uncoated), sr (blue, coated) and ms (white). Correction factors taken
% from Bonnet et al, RSI 2013. 
%
% Inputs, mandatory
%       X - signal, PSL
%       t - fading time, minutes
% Inputs, optional
%       IP type - 'tr' (default), 'ms', 'sr'
%
% Outputs
%       Xout - corrected signal, PSL
%       f - correction function
%
% Author
%   Shaun Kerr, 2016
%
%   Bonnet et al, RSI 2013
%   http://scitation.aip.org/content/aip/journal/rsi/84/10/10.1063/1.4826084
%
%%
IPtype = 'tr';
if nargin > 2
    IPtype = varargin{1};
end

switch IPtype
    
    case 'tr'
        B1 = 0.49;
        beta1 = 17.9; % min
        B2 = 0.51; 
        beta2 = 1482; % min
    case 'ms'
        B1 = 0.26;
        beta1 = 37.6; % min
        B2 = 0.74;
        beta2 = 2604; % min
    case 'sr'
        B1 = 0.49; 
        beta1 = 11.9; %min
        B2 = 0.51; 
        beta2 = 1390; % min
end

f = @(x) B1*exp(-x./beta1)+B2*exp(-x./beta2);
Xout = X/f(t);

       