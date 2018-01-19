function [Ynew, NonZeroIndex] = IonFilterCompensation(E, Y, FilterThickness, varargin)
% Returns new energy vector to account for values removed below zero. 
%   INPUTS:
%       E [MeV]
%       Y [PSL] - signal
%       FilterThickness [um]
%       (optional) 'Flag_ShowPlots'
%       (optional) 'C6'
%       
%   OUPTUTS
%       Ynew [PSL] - corrected signal
%
%   Ions & Filters
%       Protons, 50 100 250 300 um Al, 2um Mylar filters
%       Carbon 6+, 50um Al filters
%       
%   Power curve fit to downshifted ions:
%       -a*(x-b)^-d+c  
%
%   SUBFUNCTIONS:
%       psl2protons
%       psl2carbons 
%
%   Author:
%   Shaun Kerr, University of Alberta, 2018
%
%%
%---------------------------%
IonSpecies = 'protons';
Flag_ShowPlots = 1;

Flag_SavePlots = 0; 
OutputName1 = 'IPS_Filter_Step1';
OutputName2 = 'IPS_Filter_Step2';

Yrange2 = [0 1]; 
Color_Corrected = 'magenta';

% Energy at which uncertainy becomes large
% Arbitrarily selected
LowerCutoff = 63;

Erange = [0 60]; 
%---------------------------%

if ~isempty(varargin)
        if cell2mat(cellfun(@(x) strfind(x,'Flag_ShowPlots'),varargin, 'UniformOutput', 0))
            Flag_ShowPlots = 1;
        end
        if cell2mat(cellfun(@(x) strfind(x,'C6'),varargin, 'UniformOutput', 0))
            IonSpecies = 'C6';
            Erange = [0 200];
        end
end

if FilterThickness == 2
    Erange = [0 5];
end

% Make sure input is in column vectors
if isrow(E); E = E'; end
if isrow(Y); Y = Y'; end

% [Protons Carbon 6+]
switch FilterThickness
    case 2 % actually Mylar
        a = -0.07209;
        b = -0.8034;
        c = -0.01237;
    case 48.5
         a = 1.684; 
         b = 1.639;
         c = 0.1675;
         d = 0.4659; 
    case 50
         a = [1.684    -5398] ; 
         b = [1.639    -1.26];
         c = [0.1675    -4.6];
         d = 0.4659; 
    case 100
        a = 3.04; %a = 3.03;
        b = 2.74; %b = 2.74;
        c = 0.40; %c = 0.37;
        d = 0.41; %d = 0.42; 
    case 250
        % 24000 ions, 2400 ions
        a = 6.28; %a = 6.27;
        b = 5.01; %b = 5.07;
        c = 1.13;  %c = 1.22;
        d = 0.34; %d = 0.32; 
    case 300
        a = 7.24; %a = 7.35;
        b = 5.62; %b = 5.70;
        c = 1.36; %c = 1.64;
        d = 0.33; %d = 0.30; 
end

switch IonSpecies
    case 'protons'
        a = a(1);
        b = b(1);
        c = c(1); 
        if FilterThickness ~= 2
            DownshiftFunction = @(E) -a*(E-b).^(-d)+c;
        else
            DownshiftFunction = @(E) a*E.^b+c;
        end
    case 'C6'
        a = a(2);
        b = b(2);
        c = c(2); 
        DownshiftFunction = @(E) a*E.^b+c;
end
%EDownshift = flipud(Downshift(:,2));


% Calculate energy downshift 
Edownshift = DownshiftFunction(E); 
Ynew = Y; 

if Flag_ShowPlots
    %figure(1); clf; plot(E); hold on; 
    %title('1 - Energy shift'); xlabel('pixel'); ylabel('energy (MeV)')
    figure(1); clf; plot(E, Ynew); hold on
    title('1. Spectrum, Before and After Filter'); 
    xlabel('Energy (MeV)', 'FontWeight', 'bold'); 
    ylabel('PSL', 'FontWeight', 'bold')
end

% Find new proton energies
Eshifted = E + Edownshift;

if Flag_ShowPlots
    %figure(1); plot(Eshifted); 
    %legend({'original', 'shifted'})
    figure(1); plot(Eshifted, Ynew)
    xlim(Erange)
    
    Eplot = 0:0.1:200; 
    
    switch IonSpecies
        case 'protons'
            Nplot = psl2protons(Eplot);
        case 'C6'
            Nplot = psl2carbons(Eplot);
    end
    
    yyaxis right
    plot(Eplot, Nplot, ':')
    ylabel('Response curve (PSL/ion)', 'FontWeight', 'bold')
    legend({'Spectrum before filter', 'Spectrum after filter', 'Response curve'})

    if Flag_SavePlots
        SaveThesisPlot(OutputName1)
    end
end

% Remove negative energy values (protons stopped by filter)
NonZeroIndex = Eshifted>0; 
Eshifted = Eshifted(NonZeroIndex);
Ynew = Ynew(NonZeroIndex); 
NumNegative = length(E) - sum(NonZeroIndex); 

% Find IP response at these energies
switch IonSpecies
    case 'protons'
        Nout = psl2protons(Eshifted);
        Noriginal = psl2protons(E); 
    case 'C6'
        Nout = psl2carbons(Eshifted);
        Noriginal = psl2carbons(E); 
end
%size(Eshifted); size(Nout); size(Ynew)

Ynew = Ynew./Nout;
Yoriginal = Y./Noriginal; 

% Add NaN to compensate for negative values cut-off
if E(end)>E(1) % add to beginning of vector
    Ynew = [NaN(NumNegative, 1); Ynew];    
else % add to end of vector
    Ynew = [Ynew; NaN(NumNegative, 1)];
end

 if Flag_ShowPlots
     
     [~, Index2] = min(abs(E - LowerCutoff));
     
     figure(2); clf
     %plot(E, Noriginal); hold on; plot(Eshifted, Nout)
     plot(E, Yoriginal); hold on; 
     plot(E(Index2:end), Ynew(Index2:end), 'Color', Color_Corrected)
     plot(E(1:Index2), Ynew(1:Index2),':', 'Color', Color_Corrected)

     set(gca, 'YScale','linear')
     switch IonSpecies
         case 'protons'
             IonString = 'Protons'; 
         case 'C6'
             IonString = 'C^{6+}';
     end
     
     xlabel('Energy (MeV)', 'FontWeight', 'bold'); 
     ylabel([IonString ' Signal'], 'FontWeight', 'bold')
     legendstring = {'Original spectrum', 'Corrected spectrum'...
                sprintf('Corrected spectrum\n(omitted due to uncertainty)')};
     legend(legendstring)
     titlestring = sprintf('2. Results of Filter Correction, %.0f %cm Al \mu',FilterThickness,  char(181)); 
     title(titlestring)
     xlim(Erange)
     ylim(Yrange2)
     
     if Flag_SavePlots
        SaveThesisPlot(OutputName2)
    end
 end
%---------------------------------%