function Output = psl2protons(Ein, varargin)
% Use to find the absolute number of protons from a PSL signal. If only an
% energy vector is input then a corresponding response curve is output.
% If both an energy vector and PSL signal vector are given the proton
% signal is output. 
%
% The response curve is found based on one of two published results (Mancic
% et al., RSI 2008 and Bonnet et al., RSI 2013). Bonnet is the default, as 
% it is the more recent result and derived from a benchmarked GEANT4 model 
% of image plates. 
%
% Applies to TR image plates (uncoated, blue). 
%
% Usage:
%   Output = psl2protons(Ein, PSL*, paper*)
%                                               * optional
%
% Inputs, mandatory
%       Ein - energy vector or matrix (MeV)
% Inputs, optional
%       PSL - signal in PSL, vector or matrix (columns = energies)
%       paper - 'Bonnet' (default) or 'Mancic'
%
% Outputs
%       Response vector (PSL/proton) or proton signal (# protons), in 
%       ascending or descending order depending on input.
% 
% Requirements
%   Bonnet calculations require two input files in the script directory:
%       Bonnet_ExtrapolationCurve.mat 
%       Bonnet_InterpolationCurve.mat
%
% Author
%   Shaun Kerr, 2016
%
%   Mancic et al., RSI 2008
%   http://scitation.aip.org/content/aip/journal/rsi/79/7/10.1063/1.2949388
%
%   Bonnet et al, RSI 2013
%   http://scitation.aip.org/content/aip/journal/rsi/84/10/10.1063/1.4826084
%
%%
% If any input is string, take it to be paper. Otherwise it is taken as PSL
% signal
Paper = 'Bonnet';
OutputProtons = 0;
for i = 1:nargin-1
    if ischar(varargin{i})
        Paper = varargin{i};
    else
        OutputProtons = 1;
        PSL = varargin{i};
    end
end

% Make sure inputs are vertical vectors
if isrow(Ein); Ein = Ein';end
if OutputProtons
    if isrow(PSL); PSL = PSL';end
end

switch Paper
    case 'Mancic'
        % Valid for E < 2.11 MeV
        FIT1 = @(E) 0.22039*exp(-((E-1.5049).^2)/1.1842^2);

        % Valid for E = 2.11 - 20 MeV
        FIT2 = @(E) 0.33357*E.^(-0.91377);

        Nlow = FIT1(Ein(Ein<2.11));
        Nhigh = FIT2(Ein(Ein>=2.11)); 
        
    case 'Bonnet'
        % Uses previously calculated interpolated and extrapolated curves
        % to fit any energy value. 

        % Load data
        %Path = '/Users/smkerr/Dropbox/Software/Matlab/Matlab-Path-Scripts/analysis - IPS/Reference Files - rear';
        Path = fileparts(mfilename('fullpath'));
        load([Path filesep 'Bonnet_InterpolationCurve'])
        load([Path filesep 'Bonnet_ExtrapolationCurve'])
        
        Ecutoff = 40; 
        Nlow = fitresult(Ein(Ein<=Ecutoff)); 
        Nhigh = fitresult2(Ein(Ein>Ecutoff));
        
        % Load calibration for IPS
%         x = load([Path 'Bonnet_response_proton_TR_IP_40MeV.txt']);
%         E = x(:,1);
%         R = x(:,2);

         % Interpolate for 0 to 40 MeV
%         ft = 'pchipinterp';
%         [fitresult, gof] = fit( E, R, ft, 'Normalize', 'on' );
%         save([Path 'Bonnet_InterpolationCurve'], 'fitresult')
        
        % Extrapolate beyond 40 MeV
%         ft2 = fittype( 'exp2' );
%         opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
%         opts.Display = 'Off';
%         opts.StartPoint = [0.374541174828476 -0.264144333040478 0.0343676906133224 -0.00639188388493325];
%         [fitresult2, gof2] = fit( E(110:end), R(110:end), ft2, opts);
%         save([Path 'Bonnet_ExtrapolationCurve'], 'fitresult2')
end

% Flip calibration if necessary
if Ein(1) > Ein(end)
    % Energy vector is descending
    Nresponse = [Nhigh; Nlow];
else
    % Energy vector in ascending
    Nresponse = [Nlow; Nhigh];
end

if OutputProtons
    % Calculate absolute number of protons
    if size(PSL,1)>1 && size(PSL,2)>1 % PSL is matrix
        Output = PSL./Nresponse';
    else % PSL is vector
        Output = PSL./Nresponse;
    end
else
    Output = Nresponse; 
end
