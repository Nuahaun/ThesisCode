function varargout = BeamAnalysisOptions(varargin)
% BEAMANALYSISOPTIONS MATLAB code for BeamAnalysisOptions.fig
%      BEAMANALYSISOPTIONS, by itself, creates a new BEAMANALYSISOPTIONS or raises the existing
%      singleton*.
%
%      H = BEAMANALYSISOPTIONS returns the handle to a new BEAMANALYSISOPTIONS or the handle to
%      the existing singleton*.
%
%      BEAMANALYSISOPTIONS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BEAMANALYSISOPTIONS.M with the given input arguments.
%
%      BEAMANALYSISOPTIONS('Property','Value',...) creates a new BEAMANALYSISOPTIONS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before BeamAnalysisOptions_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to BeamAnalysisOptions_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help BeamAnalysisOptions

% Last Modified by GUIDE v2.5 12-Jul-2015 23:31:34

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @BeamAnalysisOptions_OpeningFcn, ...
                   'gui_OutputFcn',  @BeamAnalysisOptions_OutputFcn, ...
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


% --- Executes just before BeamAnalysisOptions is made visible.
function BeamAnalysisOptions_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to BeamAnalysisOptions (see VARARGIN)

setappdata(0  , 'hOptionsGui'    , gcf);
hOptionsGui = getappdata(0, 'hOptionsGui');       % Set handle hmain to this gui (root)

DefaultEnergy = 40; 
set(handles.edit_Energy,'string', num2str(DefaultEnergy));
setappdata(hOptionsGui,'Energy', DefaultEnergy)

DefaultPulseDuration = 1e-12;
set(handles.edit_PulseDuration,'string', num2str(1e-12*10^12));
setappdata(hOptionsGui,'PulseDuration', DefaultPulseDuration)

DefaultBkgd = 1 ;
BackgroundTypes = {'Threshold', 'RFit'};
set(handles.pop_BkgdTypes,'string', BackgroundTypes);
set(handles.pop_BkgdTypes, 'value', DefaultBkgd)
setappdata(hOptionsGui,'BackgroundType', BackgroundTypes{DefaultBkgd})

DefaultPlotType = 1 ;
PlotTypes = {'linear', 'log'};
set(handles.pop_PlotType,'string', PlotTypes);
set(handles.pop_PlotType, 'value', DefaultPlotType)
setappdata(hOptionsGui,'PlotType', PlotTypes{DefaultPlotType})

DefaultBkgdStdDev = 3; 
set(handles.edit_BkgdStdDev, 'string', num2str(DefaultBkgdStdDev))
setappdata(hOptionsGui,'BkgdStdDev', DefaultBkgdStdDev)

setappdata(hOptionsGui,'SavePlots', 0)
setappdata(hOptionsGui,'ShowPlots', 0)



DefaultMagnification = 0.343;  % H-Jet experiment, Titan 2015
set(handles.edit_Magnification ,'string', num2str(DefaultMagnification));
setappdata(hOptionsGui,'Magnification', DefaultMagnification)

% Choose default command line output for BeamAnalysisOptions
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes BeamAnalysisOptions wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = BeamAnalysisOptions_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function edit_Energy_Callback(hObject, eventdata, handles)
hOptionsGui = getappdata(0, 'hOptionsGui');       % Set handle hmain to this gui (root)
Energy = str2double(get(handles.edit_Energy, 'string'));
setappdata(hOptionsGui,'Energy',Energy);

function edit_Energy_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_PulseDuration_Callback(hObject, eventdata, handles)
hOptionsGui = getappdata(0, 'hOptionsGui');       % Set handle hmain to this gui (root)
PulseDuration = str2double(get(handles.edit_PulseDuration, 'string'))*10^-12;
setappdata(hOptionsGui,'PulseDuration',PulseDuration);

function edit_PulseDuration_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function cb_SavePlots_Callback(hObject, eventdata, handles)
hOptionsGui = getappdata(0, 'hOptionsGui');       % Set handle hmain to this gui (root)
SavePlots = get(handles.cb_SavePlots, 'value');
setappdata(hOptionsGui,'SavePlots',SavePlots);


function edit_Magnification_Callback(hObject, eventdata, handles)
hOptionsGui = getappdata(0, 'hOptionsGui');       % Set handle hmain to this gui (root)
Magnification = str2double(get(handles.edit_Magnification, 'string'));
setappdata(hOptionsGui,'Magnification',Magnification);


function edit_Magnification_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function cb_ShowPlots_Callback(hObject, eventdata, handles)
hOptionsGui = getappdata(0, 'hOptionsGui');       % Set handle hmain to this gui (root)
ShowPlots = get(handles.cb_ShowPlots, 'value');
setappdata(hOptionsGui,'ShowPlots',ShowPlots);

function pop_BkgdTypes_Callback(hObject, eventdata, handles)
hOptionsGui = getappdata(0, 'hOptionsGui');       % Set handle hmain to this gui (root)
BkgdTypes = get(handles.pop_BkgdTypes, 'string');
Bkgd = BkgdTypes{get(handles.pop_BkgdTypes, 'value')}; 
setappdata(hOptionsGui,'BkgdType',Bkgd);

function pop_BkgdTypes_CreateFcn(hObject, eventdata, handles)

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_BkgdStdDev_Callback(hObject, eventdata, handles)
hOptionsGui = getappdata(0, 'hOptionsGui');       % Set handle hmain to this gui (root)
BkgdStdDev = str2double(get(handles.edit_BkgdStdDev, 'string'));
setappdata(hOptionsGui,'BkgdStdDev',BkgdStdDev);

function edit_BkgdStdDev_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function pop_PlotType_Callback(hObject, eventdata, handles)
hOptionsGui = getappdata(0, 'hOptionsGui');      
PlotTypes = get(handles.pop_PlotType, 'string');
PlotType = PlotTypes{get(handles.pop_PlotType, 'value')}; 
setappdata(hOptionsGui,'PlotType',PlotType);

function pop_PlotType_CreateFcn(hObject, eventdata, handles)

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
