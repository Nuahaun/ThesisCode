function Output = psl2carbons(Ein, varargin)
% Use to find the absolute number of carbon ions from a PSL signal. If the  
% input is an energy vector then the output is a response curve.
% If both an energy vector and PSL signal vector are given the signal in 
% # of carbons is output. 
%
% The response curve is from Doria et al., RSI 2015. They found it to be
% independent of ion charge. 
%
% Applies to TR image plates (uncoated, blue). 
%
% Inputs, mandatory
%       Ein - energy vector or matrix (MeV)
% Inputs, optional
%       PSL - signal in PSL, vector or matrix (columns = energies)
%
% Outputs
%       Response vector (PSL/carbon) or carbon signal (# carbon), in 
%       ascending or descending order depending on input.
%
% Author
%   Shaun Kerr, 2017
%
%   Doria et al., RSI 2015
%   http://aip.scitation.org/doi/10.1063/1.4935582
%
%%
OutputIons = 0;
for i = 1:nargin-1
    OutputIons = 1;
    PSL = varargin{i};
end

% Make sure inputs are vertical vectors
if isrow(Ein); Ein = Ein';end
if OutputIons
    if isrow(PSL); PSL = PSL';end
end

%%%%%%%
% Response functions from  Doria et al., RSI 2015
% Normalized for 30 minute scan time

% Valid for E <= 73.6 MeV
FIT1 = @(E) (2.51e-3+4.56e-4*E+-8.90e-6*E.^2+4.61e-8*E.^3).*E;

% Valid for E > 73.7
FIT2 = @(E) 4.55*E.^(-0.533);

Ecutoff = 73.6;
Nlow = FIT1(Ein(Ein<=Ecutoff));
Nhigh = FIT2(Ein(Ein>Ecutoff)); 
        
% Flip calibration if necessary
if Ein(1) > Ein(end)
    % Energy vector is descending
    Nresponse = [Nhigh; Nlow];
else
    % Energy vector in ascending
    Nresponse = [Nlow; Nhigh];
end

if OutputIons
    % Calculate absolute number of protons
    if size(PSL,1)>1 && size(PSL,2)>1 % PSL is matrix
        Output = PSL./Nresponse';
    else % PSL is vector
        Output = PSL./Nresponse;
    end
else
    Output = Nresponse; 
end
