%%
%
%   Shaun Kerr, 2017
%   smkerr@gmail.com
%
%% Load Dispersion
% DispersionPath = '/Users/smkerr/Dropbox/Software/Matlab/Matlab-Path-Scripts/analysis - IPS/Reference Files - rear/Rear IP dispersion IPS1.xlsx';
% [~, ~, raw] = xlsread(DispersionPath,'Sheet1');
% raw = raw(66:end,:);
% temp = reshape([raw{:}],size(raw));
% Dispersion = [temp(:,1), temp(:,3)];
% Energies = Dispersion(:,2);

%% Create energy vector and save as TRIM.dat
disp('---------')
Filename = '/Users/smkerr/Desktop/SRIM Files/TRIM.dat';
%Filename = '/Users/smkerr/Desktop/Test_IPS_Data/TRIM.dat';

Emin = 0.1; 
Emax = 3; 
Estep = 0.05; 
AtomicNumber = 1; 
RepeatNum = 150; 
Description = 'Proton energy spectra for downshift in aluminized mylar ';

Energies = Emin:Estep:Emax; Energies = Energies';

Energies = repmat(Energies, RepeatNum,1);
Energies = Energies*1e6; % convert to MeV
NumPart = length(Energies);
fprintf('Number of particles = %.0f', NumPart)
AtomicNum = ones(NumPart,1)*AtomicNumber;
Depth = zeros(NumPart,1);
Y = zeros(NumPart,1);
Z = zeros(NumPart,1);
CosX = ones(NumPart,1);
CosY = zeros(NumPart,1);
CosZ = zeros(NumPart,1);

Data = [AtomicNum,Energies, Depth, Y, Z, CosX, CosY, CosZ];

%Data = [1 1e6       0 0 0 1 0 0;...
%        2 1.564e4   0 0 0 1 0 0;...
%        3 1.5e6     0 0 0 1 0 0];


%
SaveTRIMFile(Filename, Data, Description);
disp('Write done.')