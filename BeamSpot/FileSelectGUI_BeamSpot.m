function varargout = FileSelectGUI_BeamSpot(varargin)
% FILESELECTGUI_BEAMSPOT MATLAB code for FileSelectGUI_BeamSpot.fig
%      FILESELECTGUI_BEAMSPOT, by itself, creates a new FILESELECTGUI_BEAMSPOT or raises the existing
%      singleton*.
%
%      H = FILESELECTGUI_BEAMSPOT returns the handle to a new FILESELECTGUI_BEAMSPOT or the handle to
%      the existing singleton*.
%
%      FILESELECTGUI_BEAMSPOT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FILESELECTGUI_BEAMSPOT.M with the given input arguments.
%
%      FILESELECTGUI_BEAMSPOT('Property','Value',...) creates a new FILESELECTGUI_BEAMSPOT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before FileSelectGUI_BeamSpot_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to FileSelectGUI_BeamSpot_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help FileSelectGUI_BeamSpot

% Last Modified by GUIDE v2.5 19-Aug-2015 21:53:29

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @FileSelectGUI_BeamSpot_OpeningFcn, ...
                   'gui_OutputFcn',  @FileSelectGUI_BeamSpot_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before FileSelectGUI_BeamSpot is made visible.
function FileSelectGUI_BeamSpot_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to FileSelectGUI_BeamSpot (see VARARGIN)

setappdata(0  , 'hMainGui'    , gcf);
hMainGui = getappdata(0, 'hMainGui');       % Set handle hmain to this gui (root)

Paths.Data = '/Users/smkerr/Dropbox/Grad School/Experiments/LLNL/16 - Summer 2015/Data/BeamSpot/Shot18withoutJet/'; 
%DefaultPath = '/Users/smkerr/Dropbox/Grad School/Experiments/LLNL/16 - Summer 2015/Data/BeamSpot/Shots/';
%Paths.Output = '/Users/smkerr/Dropbox/Grad School/Experiments/LLNL/16 - Summer 2015/Analysis/Beam Spots/';
Paths.Output = Paths.Data; 
setappdata(hMainGui,'Paths',Paths);

GetFiles(handles, 'csv')
setappdata(hMainGui,'CurrentRow',1);

% Launch Options GUI
BeamAnalysisOptions


% Choose default command line output for FileSelectGUI_BeamSpot
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

set(handles.table_Files,'CellSelectionCallback',{@TableCallBack, handles})


% UIWAIT makes FileSelectGUI_BeamSpot wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = FileSelectGUI_BeamSpot_OutputFcn(hObject, eventdata, handles) 
% Get default command line output from handles structure
varargout{1} = handles.output;


function GetFiles(handles, identifier)
hMainGui = getappdata(0,'hMainGui');
Paths = getappdata(hMainGui, 'Paths');
Files = dir([Paths.Data '*' identifier]);
FileArray = {Files.name}';
set(handles.table_Files, 'data', FileArray)
NumFiles = size(FileArray, 1);
set(handles.edit_NumFiles, 'string', num2str(NumFiles)) 
setappdata(hMainGui,'FileArray',FileArray);

function btn_InputPath_Callback(hObject, eventdata, handles)
hMainGui = getappdata(0, 'hMainGui'); 
Paths = getappdata(hMainGui, 'Paths'); 
Paths.Data = GetDirectory;
setappdata(hMainGui, 'Paths', Paths)
FitString = 'csv';
GetFiles(handles, FitString)


function btn_AnalyzeFile_Callback(hObject, eventdata, handles)
hOptionsGui = getappdata(0, 'hOptionsGui');       % Set handle hmain to this gui (root)
hMainGui = getappdata(0, 'hMainGui');       % Set handle hmain to this gui (root)

Paths = getappdata(hMainGui,'Paths');
Row = getappdata(hMainGui,'CurrentRow');
FileArray = get(handles.table_Files, 'data');
FileName = FileArray{Row}; 
Data = load([Paths.Data FileName]);

% Temp line
%Data = Data.A; 

[IntensityDist, EncircledEnergy] = AnalyzeFile(hMainGui, hOptionsGui, Data, FileName);

    

function [IntensityDist, EncircledEnergy] = AnalyzeFile(hMainGui, hOptionsGui, Data, FileName)
Paths = getappdata(hMainGui,'Paths'); 


magnification = getappdata(hOptionsGui, 'Magnification'); 
Energy = getappdata(hOptionsGui, 'Energy'); 
PulseDuration = getappdata(hOptionsGui, 'PulseDuration'); 
ShowPlots = getappdata(hOptionsGui, 'ShowPlots');
BkgdType = getappdata(hOptionsGui, 'BkgdType'); 
BkgdStdNum = getappdata(hOptionsGui, 'BkgdStdDev');
SavePlots = getappdata(hOptionsGui, 'SavePlots');
PlotType = getappdata(hOptionsGui, 'PlotType');

if SavePlots
[IntensityDist, EncircledEnergy]= BeamSpotFunction(Data,...
        magnification,...
        'FileName', FileName,...
        'Energy',Energy, ...
        'pulseDuration', PulseDuration, ...
        'BkgdType', BkgdType,...
        'BkgdStdNum', BkgdStdNum,...
        'ShowPlots', ShowPlots,...
        'PlotType', PlotType,...
        'Paths.Output', Paths.Output);
else
    [IntensityDist, EncircledEnergy]= BeamSpotFunction(Data,...
        magnification,...
        'FileName', FileName,...
        'Energy',Energy, ...
        'pulseDuration', PulseDuration, ...
        'BkgdType', BkgdType,...
        'BkgdStdNum', BkgdStdNum,...
        'PlotType', PlotType,...
        'ShowPlots', ShowPlots);
end

function TableCallBack(hObj,evt, handles)
%evt.Indices
CurrentRow =  evt.Indices(:,1); 
hMainGui = getappdata(0, 'hMainGui');       % Set handle hmain to this gui (root)
setappdata(hMainGui,'CurrentRow',CurrentRow);
%display(['Table selected - row ' num2str(Curr%entRow)])

function UpdateTable(handles)
hMainGui = getappdata(0, 'hMainGui');       % Set handle hmain to this gui (root)
FileArray = getappdata(hMainGui,'FileArray');
RefString = get(handles.edit_ParseString, 'string');

% No string given, show all files
if isempty(RefString)
    set(handles.table_Files, 'data', FileArray)
    NumFiles = size(FileArray,1)
    set(handles.edit_NumFiles, 'string', num2str(NumFiles))
    return
end

Index = strfind(FileArray, RefString);
Index = find(~cellfun(@isempty, Index));

NewFileArray = FileArray(Index); 
NumFiles = size(NewFileArray,1);

set(handles.table_Files, 'data', NewFileArray)
set(handles.edit_NumFiles, 'string', num2str(NumFiles))


function edit_ParseString_Callback(hObject, eventdata, handles)
UpdateTable(handles)

function edit_ParseString_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_NumFiles_Callback(hObject, eventdata, handles)


function edit_NumFiles_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function btn_AnalyzeAll_Callback(hObject, eventdata, handles)


function AverageCalculations(handles)
hOptionsGui = getappdata(0, 'hOptionsGui');       % Set handle hmain to this gui (root)
hMainGui = getappdata(0, 'hMainGui');       % Set handle hmain to this gui (root)

Paths = getappdata(hMainGui,'Paths');
Row = getappdata(hMainGui,'CurrentRow');
FileArray = get(handles.table_Files, 'data');
BkgdType = getappdata(hOptionsGui, 'BkgdType'); 

Fig20 = figure(20); clf
Fig21 = figure(21); clf

legendstring = {}; 
for i = 1:size(Row,1)
    FileName = FileArray{Row(i)}
    Data = load([Paths.Data FileName]);
    [IntensityDist{i}, EncircledEnergy{i}] = AnalyzeFile(hMainGui, hOptionsGui, Data, FileName(i));
    figure(Fig20)
    PlotEncircledEnergy(EncircledEnergy{i}(:,1),EncircledEnergy{i}(:,2),1); hold on 
    
    figure(Fig21)
    PlotIntensityDistribution(IntensityDist{i}(:,1),IntensityDist{i}(:,2),1); hold on 

    legendstring = [legendstring, {FileName}];
end
legend(legendstring)

tempEE = [EncircledEnergy{:}]; 
Ravg = mean(tempEE(:,[1:2:end]),2); 
Eavg = mean(tempEE(:,[2:2:end]),2); 

Emin = min(tempEE(:,[2:2:end]),[],2); 
Emax = max(tempEE(:,[2:2:end]),[],2); 

tempID = [IntensityDist{:}]; 
EfracAvg = mean(tempID(:,[1:2:end]),2); 
IntAvg = mean(tempID(:,[2:2:end]),2); 

IntMin = min(tempID(:,[2:2:end]),[],2); 
IntMax = max(tempID(:,[2:2:end]),[],2); 

% Plot encircled energy
Fig22 = figure(22); clf
PlotEncircledEnergy(Ravg, Eavg)
title1 = get(gca, 'title');

a = strfind(FileName, 'sh');
OutputName = FileName(a(1):a(1)+3); 

titlestring = title1.String; 
titlestring = [OutputName ' - ' titlestring ' (' BkgdType ')']; 
title(titlestring)

patch([Ravg; flipud(Ravg)], [Emax; flipud(Emin)],...
    'blue', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); 

% Plot intensity distribution
Fig23 = figure(23); clf
PlotIntensityDistribution(EfracAvg, IntAvg)
title2 = get(gca, 'title');
xlim([0 1])

titlestring2 = title2.String; 
titlestring2 = [OutputName ' - ' titlestring2 ' (' BkgdType ')']; 
title(titlestring2)

ZeroIndex = min(find(IntMax == 0,1,'first'),find(IntMin == 0,1,'first'))-1;
patch([EfracAvg(1:ZeroIndex); flipud(EfracAvg(1:ZeroIndex))], [IntMax(1:ZeroIndex); flipud(IntMin(1:ZeroIndex))],...
    'blue', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); 


if get(handles.cb_SavePlots, 'value')
    Paths = getappdata(hMainGui, 'Paths'); 

    % EPS extension gives a nice vector image but it's huge (1.8 MB)
    %ext = 'eps'; % with Open GLp
    ext = 'png'; 
    fprintf('\n%s%s%s%s%s', 'Saving ', OutputName,'.', ext, '...  ')
    export_fig(Fig22, [Paths.Output OutputName '_' BkgdType '_AverageEnclosedEnergy'], ['-' ext], '-opengl')
    export_fig(Fig23, [Paths.Output OutputName '_' BkgdType '_AverageEIntensityDistribution'], ['-' ext], '-opengl')

    %export_fig(Fig10, [Paths.Output OutputName], '-eps', '-painters', '-nocrop')
    % '-native'
    %print(Fig10, '-depsc','-tiff','-painters', '-r300', [Paths.Output OutputName])
    savevar = [IntensityDist, EncircledEnergy]; 
    save([Paths.Output OutputName '_' BkgdType '_Average.mat'], 'savevar')
    fprintf('%s\n', 'Saved!')
end


function btn_Average_Callback(hObject, eventdata, handles)
clc
% Disable interface during calculation
if strcmp(get(handles.btn_Average,'String'),'Running')
    % Stopping simulation
    switchButtons(handles,'on')
    set(handles.btn_Average,'String','Average')
    return
else
    display('off!')
    set(handles.btn_Average,'String','Running')
    %set(findobj('Enable','on'),'Enable','off')
    switchButtons(handles,'off')
    %set(handles.btn_Average,'Enable','on')
    drawnow()
end

try
    AverageCalculations(handles)
catch
    warning('Problem with a spot image.');

end
% Enable interface after calculation
if strcmp(get(handles.btn_Average,'String'),'Running')
    % Stopping simulation
    switchButtons(handles,'on')
    set(handles.btn_Average,'String','Average')
    return
else
    set(handles.btn_Average,'String','Running')
    %set(findobj('Enable','on'),'Enable','off')
    switchButtons(handles,'off')
end

function switchButtons(hs,onoff)
set([hs.btn_AnalyzeFile hs.btn_AnalyzeAll hs.btn_Average],...
    'Enable',onoff)

function btn_OutputPath_Callback(hObject, eventdata, handles)


function cb_SavePlots_Callback(hObject, eventdata, handles)
