function [IntensityDist, EncircledEnergy] =  BeamSpotFunction(Spot, Mag, varargin)
% Analyzes laser spot,  giving encircled energy vs. radius, R50 and R80
% plot and intensity distribution. Two background subtraction methods are
% possible:
%       1) Thresholding, with all pixels < background mean + X standard
%       deviations set to zero
%       2) Fit to encircled energy, and subtraction
%
% Fit energy options are peak or center of mass
%
% Requires PlotEncircledEnergy function
%
%   [IntensityDist, ...
%           Column 1 = cumulative fraction of incident energy
%           Column 2 = intensity (W/cm2)
%   EncircledEnergy] 
%           Column 1 = radius (um)
%           Column 2 = relative encircled energy
%      = BeamSpotFunction(
%    REQUIRED:
%          Spot,.............. 2D matrix of spot iamge
%          Magnification,..... in um per pixel
%    OPTIONAL:
%          'energy',.......... [Energy] (in J, default is 40 J) 
%          'pulseduration', .. [PulseDuration] (in sec, default is 1ps)
%          'outputpath', ..... OutputPath, and save plots to here
%          'showplots', ...... [0 or 1 (default)] 
%          'plotstyle'         {'linear', 'log'}
%          'filename',........ FileName (for title of plot, output file)
%          'bkgdtype',........ {'RFit', 'Threshold'}
%          'bkgdstdnum',...... [number] (number of standard deviations
%          above background used for threshold)
%
% v.2
%  - inputParser for input
%
%   Shaun Kerr, 2017
%   smkerr@gmail.com.ca
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parse input variables, set defaults
set(0,'defaultaxesfontsize', 16)
p = inputParser;

defaultEnergy = 50; 
defaultPulseDuration = 1*10^-12; % ps, again nominal value

defaultOutputPath = {};
defaultFileName = '';

defaultBkgdType = 'Threshold'; 
validBkgdType = {'Threshold', 'RFit'};
checkBkgdType = @(x) any(validatestring(x,validBkgdType));

defaultPlotType = 'linear'; 
validPlotType = {'linear', 'log'};
checkPlotType = @(x) any(validatestring(x,validPlotType));


defaultShowPlots = 1; 
defaultBkgdStdNum = 3; 

addRequired(p, 'Spot', @isnumeric);
addRequired(p, 'Mag', @isnumeric);
addOptional(p,'OutputPath',defaultOutputPath);
addOptional(p,'FileName',defaultFileName,@ischar);
addOptional(p,'ShowPlots',defaultShowPlots,@isnumeric);

addOptional(p,'PlotType',defaultPlotType,checkPlotType);
addOptional(p,'BkgdType',defaultBkgdType,checkBkgdType);
addOptional(p,'BkgdStdNum',defaultBkgdStdNum,@isnumeric);

addParameter(p,'Energy',defaultEnergy,@isnumeric);
addParameter(p,'PulseDuration',defaultPulseDuration,@isnumeric);

parse(p,Spot,Mag,varargin{:});

SubplotRows = 2;
SubplotColumns = 4;
Spot1 = double(Spot); 
umPerPixel = Mag;

%MeasurementType = 'Radius';
MeasurementType = 'Diameter';

% Real magnification: 0.343
%  0.7154


% What are these for? Number of pixels?
% not sure why it was hard-coded before
xPixels = size(Spot1, 2)
yPixels = size(Spot1, 1)

%xPixels = 892;
%yPixels = 671; 

%yPixels = 480;
%xPixels = 640; 

%% Display raw image 
if p.Results.ShowPlots
    Fig1 = figure(1); clf; 
    imagesc(Spot1)
    colorbar
    axis equal
    title('Raw Image')
    set(gcf,'Color','White'); % White background
    ylim([0 yPixels])
end

%% Find background from whole image
margin = 75; 
%slice.top = Spot1(1:margin,:);
%slice.bottom = Spot1(end-margin:end,:);
slice.top = Spot1(:,1:margin);
slice.bottom=Spot1(:,end-margin:end); 


mean1 = mean(mean(slice.top));
mean2 = mean(mean(slice.bottom));
bkgd = mean([mean1 mean2]); 
StdBkgd = mean([mean(std(slice.bottom)), mean(std(slice.top))]);

%fprintf('%s\t%f\n', 'Background = ', bkgd)
%fprintf('%s\t%f\n', 'Std Dev = ', StdBkgd)

if p.Results.ShowPlots
    Fig2 = figure(2); clf;
    imagesc(Spot1)
    colormap(gca,'parula')
    title('Background Selection')
    patch([0,0, margin, margin],[0 yPixels yPixels 0],'r', 'FaceAlpha', 0.45)
    patch([xPixels-margin,xPixels-margin, xPixels, xPixels],[0 yPixels yPixels 0],'r', 'FaceAlpha', 0.45)
    %patch([0,0, 640, 640],[0 margin margin 0],'r', 'FaceAlpha', 0.45)
    %patch([0,0, 640, 640],[480 480-margin 480-margin 480],'r', 'FaceAlpha', 0.45)
    set(gcf,'Color','White'); % White background
end

%%  Crop around max pixels
% Find X and Y location of max pixel
Spot2 = Spot1; 
[~, Indpk] = max(Spot2(:));
[Ypk, Xpk]= ind2sub([size(Spot1,1) size(Spot1,2)], Indpk);
%Xpk = 238;
%Ypk = 186; 


% Crop image
BkgdStdNum = p.Results.BkgdStdNum;
switch p.Results.BkgdType
    case 'RFit'
        CropPixels = 100; 
        %CropPixels = 180; 
    case 'Threshold'
        CropPixels = 100; % Normal crop size
        %CropPixels = 50; 
        Spot2(Spot2<bkgd+BkgdStdNum*StdBkgd) = 0; 
end

% Check that spot is away from edges so crop can occur
if (Ypk-CropPixels)<0 
    error('Spot is too close to left side of image!')
    return   
elseif (Xpk-CropPixels)<0 
   error('Spot is too close to top of image!')
   return    
elseif (Ypk+CropPixels) > yPixels 
    error('Spot is too close to right side of image!')
   return  
elseif (Xpk+CropPixels) > xPixels
   error('Spot is too close to bottom of image!')
   return
end
    
Spot3 = Spot2(Ypk-CropPixels:Ypk+CropPixels,Xpk-CropPixels:Xpk+CropPixels); 

if strcmp(p.Results.BkgdType, 'Threshold')
    Fig13 = figure(13); clf; 
    Img3 = imagesc(Spot3);
    %colormap(gca,'parula')
    axis square
    colorbar
    colormap jet
    ylabel('pixels')
    xlabel('pixels')

    title(['Cropped, Background Threshold = Mean + ' num2str(BkgdStdNum) '\sigma'])
    set(gcf,'Color','White'); % White background
    %colormapeditor
    mycmap = colormap(gca);
    %save('MyColormaps','mycmap')
    mycmap(1,:) = [1 1 1];
    colormap(gca,mycmap);
end

%load('MyColormaps','mycmap')

%colormap(mycmap)
%set(gcf,'colormap',mycmap)
%drawnow
%mycmap = colormap(ax)
%d3 = Spot3;
%d3(Spot3<0) = 5/0
%imagesc(d3)
%colorbar

%% Find enclosed energy in circles of arbitrary radius
%PeakType = 'centerofmass';
PeakType = 'peak'; 
% Find X and Y location of max pixel
[~, Indpk] = max(Spot3(:));
[Ypk, Xpk]= ind2sub([size(Spot3,1) size(Spot3,2)], Indpk);
%Xpk = 238;
%Ypk = 186; 
switch PeakType
    case 'centerofmass'
        % Find center of mass
        TempSpot = Spot3;
        TempSpot(TempSpot<bkgd+5*StdBkgd) = 0; 
        if p.Results.ShowPlots
            figure(11)
            imagesc(TempSpot)
            title('Center of Mass')
        end
        A = centerOfMass(TempSpot);
        XcenterPixel = A(2);
        YcenterPixel = A(1); 
    case 'peak'
        XcenterPixel = Xpk;
        YcenterPixel = Ypk; 
end

%XcenterPixel = 227;
%YcenterPixel = 150;
XCenterUm = XcenterPixel*umPerPixel;
YCenterUm = YcenterPixel*umPerPixel;

XSize = size(Spot3,1); 
YSize = size(Spot3,2);
[X,Y] = meshgrid(1:1:XSize, 1:1:YSize);
x = [1:XSize]*umPerPixel;
y = [1:YSize]*umPerPixel;

TotalSig = sum(sum(Spot3)); 
    Ind = 1;
R = 0:0.1:50;
for i = R
    dist = sqrt((X-XcenterPixel).^2 + (Y-YcenterPixel).^2)*umPerPixel;
    EnclosedFrac(Ind) = sum(Spot3(find(dist<i)))/TotalSig;
    Ind = Ind+1;
end

%% Fit to r^2, subtract fit
switch p.Results.BkgdType
    case 'RFit'
        i1 = round(length(R)*0.5); 
        [pFit,S] = polyfit(R(i1:end).^2,EnclosedFrac(i1:end), 1);
        EnclosedEFit = polyval(pFit, R.^2);

       if p.Results.ShowPlots
            Fig4 = figure(4); clf; 
            plot(R.^2*pi, EnclosedFrac); hold on
            plot(R(i1:end).^2*pi, EnclosedFrac(i1:end), 'go')
            plot(R.^2*pi, EnclosedEFit, 'r--', 'LineWidth', 2)
            grid on
            xlabel('Circle Area (um^2)')
            ylabel('Encircled Relative Energy')
            title('Encircled Energy')
            set(gcf,'Color','White'); % White background
            legend({'Encirlced Energy','Fit Range','Fit'})
            text(max(R.^2)/2, double(max(EnclosedEFit)/2), [num2str(pFit(1)) 'x + ' num2str(pFit(2))],...
                    'Color', 'black','FontSize', 18)
        end
        CorrectedEFrac = EnclosedFrac-EnclosedEFit;
        CorrectedEFrac = CorrectedEFrac+abs(min(CorrectedEFrac));
        CorrectedEFrac = CorrectedEFrac/max(CorrectedEFrac); 
    case 'Threshold'
        CorrectedEFrac = EnclosedFrac;
end

EncircledEnergy = [R; CorrectedEFrac]';



%% Plot rings of enclosed energy
switch p.Results.BkgdType
    case 'RFit'
        CircleColor = 'red';
        TextColor = 'red';
        BkgdString = 'Bkgd = Enclosed Energy Fit'; 
    case 'Threshold'
        CircleColor = 'red';
        TextColor = 'red';
        %TextColor = [0.851 0.325 0.098];
        BkgdString = ['Bkgd Threshold = Mean + ' num2str(BkgdStdNum) 'sigma']; 

end
        
% Plot spot with overlayed circles
Fig10 = figure(10); clf; 
h(2) = subplot(SubplotRows,SubplotColumns,[1 2 5 6]);

switch p.Results.PlotType
    case 'log'
        imagesc(x,y,log(Spot3/max(max(Spot3)))); hold on
        titlestring = 'Log'; 
        colormaptype = 'parula'; 
    case 'linear' 
        imagesc(x,y,Spot3/max(max(Spot3))); hold on
        titlestring = 'Linear'; 
        colormaptype = 'jet'; 
end
plot(XCenterUm, YCenterUm, 'wo');
%colormap(gca,mycmap)
set(gcf,'Color','White'); % White background
eval(['colormap ' colormaptype])
axis square
xlabel('x (um)')
ylabel('y (um)')

switch MeasurementType
    case 'Radius'
        titlestring2 = ' plot with 50% & 80% enclosed energy radii';
        nstring = 'R'; 
    case 'Diameter'
        titlestring2 = ' plot with 50% & 80% enclosed energy diameters';
        nstring = 'D'; 
end

titlestring = sprintf('%s', [titlestring titlestring2]);
%titlestring = sprintf('%s\n%s', ' R50 and R80 spots', FileName);
title(titlestring, 'Interpreter', 'none')
colorbar
LimUM = 20; 
xlim([XCenterUm-LimUM XCenterUm+LimUM])
ylim([YCenterUm-LimUM YCenterUm+LimUM])

Ind50 = find(CorrectedEFrac>0.5, 1, 'first');
if ~isempty(Ind50)
    R50 = R(Ind50); 
    switch MeasurementType
        case 'Radius'
            num50 = R50; 
        case 'Diameter' 
            num50 = R50*2; 
    end
    viscircles([XCenterUm YCenterUm], R50,'EdgeColor',CircleColor, 'LineStyle', '--', 'LineWidth', 4);
    text(XCenterUm-4*umPerPixel, YCenterUm+R50+4*umPerPixel, [num2str(num50) 'um'],...
        'Color', TextColor,'FontSize', 22, 'FontWeight', 'bold','BackgroundColor', 'white','EdgeColor', 'black')

    
    
    Ind80 = find(CorrectedEFrac>0.8, 1, 'first');
    
    r50string = [nstring '50 = ' num2str(num50) ' um'];
    boxstring = r50string;
     if ~isempty(Ind80)
         R80 = R(Ind80); 
        switch MeasurementType
            case 'Radius'
                num80 = R80; 
            case 'Diameter' 
                num80 = R80*2; 
        end
         viscircles([XCenterUm YCenterUm], R80,'EdgeColor',CircleColor, 'LineStyle', '--', 'LineWidth', 4);

        text(XCenterUm-4*umPerPixel, YCenterUm+R80+4*umPerPixel, [num2str(num80) 'um'],...
            'Color', TextColor,'FontSize', 22, 'FontWeight', 'bold','BackgroundColor', 'white','EdgeColor', 'black')

        r80string = [nstring '80 = ' num2str(num80) ' um'];
        boxstring = sprintf('%s\n%s', boxstring, r80string);
     end

     %text(XCenterUm+7, YCenterUm-17, boxstring,...
     %   'Color', 'red','FontSize', 18,'BackgroundColor', 'white','EdgeColor', 'black')
end

if strcmp(p.Results.BkgdType, 'Threshold')
        colormap(gca,mycmap);
end

% Modify position to fit
h2Pos = get(h(2), 'Position');
set(h(2), 'position', [h2Pos(1), h2Pos(2), h2Pos(3)*0.85, h2Pos(4)] );

%% Plot enclosed energy fraction on main plot
figure(Fig10); 
set(Fig10, 'Position', [175 222 1452 726]);
h(1) = subplot(SubplotRows,SubplotColumns,[3 4]);

PlotEncircledEnergy(R, CorrectedEFrac, MeasurementType)

%{
plot(R, CorrectedEFrac, 'LineWidth', 2); 
grid on
xlim([0 30])
xlabel('Radius (um)')
ylabel('Relative Encircled Energy')
title('Encircled Energy vs. Radius')
set(gcf,'Color','White'); % White background

% Draw lines at R50 and R80
line([R50 R50 0], [0 CorrectedEFrac(Ind50) CorrectedEFrac(Ind50)],'Color', 'red', 'LineWidth', 2,'LineStyle', '--' )
line([R80 R80 0], [0 CorrectedEFrac(Ind80) CorrectedEFrac(Ind80)],'Color', 'red', 'LineWidth', 2,'LineStyle', '--' )
text(1, CorrectedEFrac(Ind50)+0.05, [num2str(R50) ' um'],...
    'Color', 'red','FontSize', 18)
text(1, CorrectedEFrac(Ind80)+0.05, [num2str(R80) ' um'],...
    'Color', 'red','FontSize', 18)
%}

%% Calculate intensity distribution based on circles


%% Find spot to calculate intensity in
Ind99 = find(CorrectedEFrac>0.99, 1, 'first');
R99 = R(Ind99); 
IntDistSpot = Spot3; 
IntDistSpot = IntDistSpot - bkgd; 

IntDistSpot(dist>R99)=0;
IntDistSpot(IntDistSpot<0) = 0; 

if p.Results.ShowPlots
    figure(9)
    imagesc(log(IntDistSpot))
    colorbar
    colormap jet
    title('Log plot of spot enclosing 99% energy')
    axis equal
end
%% Calculate intensity
Energy = p.Results.Energy;
PulseDuration = p.Results.PulseDuration; 
% Sort pixels in descending order
OrderedCounts = sort(IntDistSpot(:),'descend');

 % Remove negative counts
[I1, b2]  = find(OrderedCounts<0, 1, 'first');
if isempty(I1) % no negative counts
    PixelIndex = 1:length(OrderedCounts); 
else
    PixelIndex = 1:I1-1;
end

PosCounts = OrderedCounts(PixelIndex);
 
Counts = PosCounts;

CumCounts = cumsum(Counts);
Sum1 = sum(Counts); 
FracBeamSpot = CumCounts/Sum1; 

%umPerPixel = 0.35; 
cmPerPixel = umPerPixel/10^4; 
%areaPerPixel_um = umPerPixel^2; 
areaPerPixel_cm = cmPerPixel^2; 
%PixelArea = PixelIndex*areaPerPixel;
%PixelRadiusUm = sqrt(PixelArea/pi);

A1 = Counts/Sum1;
IntDist = A1*Energy/PulseDuration/areaPerPixel_cm;
IntensityDist = [FracBeamSpot, IntDist]; 

Imin = 10^17; 
Imax = 2*10^20; 
ystring = 'Intensity (W/cm^{2})';
Ymin = Imin; 
Ymax = Imax; 
Icounts = round(length(Counts)/4); 

% Intensity distribution plot on main plot
figure(Fig10); 
subplot(SubplotRows,SubplotColumns,[7 8])
suptitle([p.Results.FileName ' (' BkgdString ')'])
PlotIntensityDistribution(FracBeamSpot(1:Icounts),IntDist(1:Icounts), 0, Energy, PulseDuration)
%{
semilogy(FracBeamSpot(1:Icounts), IntDist(1:Icounts), 'LineWidth', 2)
set(gcf,'Color','White'); % White background
xlabel('Cumulative Fraction of Incident Energy')
ylabel(ystring)
grid on; hold on;
ylim([Ymin Ymax])
titlestring2 = sprintf('%s%0.2f%s%0.1f%s','Pixel Intensity Distribution (', Energy, ' J, ',PulseDuration*10^12, ' ps)');
title(titlestring2)
tempstring = sprintf('%s%0.1i', 'Imax = ', max(IntDist(1:Icounts)));
text(0.65, max(IntDist(1:Icounts)), tempstring,...
    'Color', 'black','FontSize', 18,...
    'BackgroundColor', 'white','EdgeColor', 'black')
IndexIntDist50 = find(FracBeamSpot<0.5, 1, 'last');
IntDist50 = IntDist(IndexIntDist50);
line([0.5 0.5 0], [Imin IntDist50 IntDist50],'Color', 'red', 'LineWidth', 2,'LineStyle', '--' )
I50string = sprintf('%0.1e',IntDist50);
text(0.01, IntDist50*0.5, I50string,...
    'Color', 'red','FontSize', 18)
%}
%FitRange = Icounts-35000:Icounts; 
%semilogy(FracBeamSpot(FitRange), IntDist(FitRange), 'r-', 'LineWidth', 4)
%[p1, S1] = polyfit(FracBeamSpot(FitRange), log(IntDist(FitRange)),1);
%hold on;
%semilogy(FracBeamSpot, 10.^polyval(p1, FracBeamSpot))

LaserWavelength = 527; %nm
ystring = ['a_0 for \lambda = ' num2str(LaserWavelength) ' nm'];
Ymin = 0.8549*sqrt( (LaserWavelength/1000)^2*(Imin/1e18) );
Ymax = 0.8549*sqrt( (LaserWavelength/1000)^2*(Imax/1e18) );
a0 = 0.8549*sqrt( (LaserWavelength/1000)^2*(IntDist./1e18) );

% Plot a0
if p.Results.ShowPlots
    Fig8 = figure(8); clf; 
    semilogy(FracBeamSpot, a0)
    set(gcf,'Color','White'); % White background
    xlabel('Cumulative Fraction of Incident Energy')
    ylabel(ystring)
    grid on
    ylim([Ymin Ymax])
    title('Pixel Intensity Distribution')
    %plot(FracBeamSpot, PixelRadiusUm)
end

if ~isempty(p.Results.OutputPath)
    a = strfind(p.Results.FileName, '.');
    OutputName = p.Results.FileName(1:a(1)-1); 
    % EPS extension gives a nice vector image but it's huge (1.8 MB)
    %ext = 'eps'; % with Open GLp
    ext = 'png'; 
    fprintf('\n%s%s%s%s%s', 'Saving ', OutputName,'.', ext, '...  ')
    export_fig(Fig10, [p.Results.OutputPath OutputName '_AnalyzedSpot'], ['-' ext], '-opengl')
    %export_fig(Fig10, [OutputPath OutputName], '-eps', '-painters', '-nocrop')
    % '-native'
    %print(Fig10, '-depsc','-tiff','-painters', '-r300', [OutputPath OutputName])
    fprintf('%s\n', 'Saved!')
end