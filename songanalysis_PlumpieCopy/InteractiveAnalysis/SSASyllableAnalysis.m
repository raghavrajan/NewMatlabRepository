function varargout = SSASyllableAnalysis(varargin)
% SSASYLLABLEANALYSIS MATLAB code for SSASyllableAnalysis.fig
%      SSASYLLABLEANALYSIS, by itself, creates a new SSASYLLABLEANALYSIS or raises the existing
%      singleton*.
%
%      H = SSASYLLABLEANALYSIS returns the handle to a new SSASYLLABLEANALYSIS or the handle to
%      the existing singleton*.
%
%      SSASYLLABLEANALYSIS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SSASYLLABLEANALYSIS.M with the given input arguments.
%
%      SSASYLLABLEANALYSIS('Property','Value',...) creates a new SSASYLLABLEANALYSIS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SSASyllableAnalysis_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SSASyllableAnalysis_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SSASyllableAnalysis

% Last Modified by GUIDE v2.5 18-Feb-2013 23:17:36

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SSASyllableAnalysis_OpeningFcn, ...
                   'gui_OutputFcn',  @SSASyllableAnalysis_OutputFcn, ...
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


% --- Executes just before SSASyllableAnalysis is made visible.
function SSASyllableAnalysis_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SSASyllableAnalysis (see VARARGIN)

% Choose default command line output for SSASyllableAnalysis
handles.output = hObject;

% Initialise some variables using either default values or values that have
% been passed while calling the function (RR 07 Jan 2013)

if (nargin >= 1)
    Temp = varargin{1};
    handles.SSASA = Temp.SSA;
    handles.SSASA.PreTime = 0.2;
    handles.SSASA.PostTime = 0.2;
    TempLabels = [];
    if (isfield(handles.SSASA.UnDirFileInfo, 'NoteLabels'))
        for i = 1:length(handles.SSASA.UnDirFileInfo.NoteLabels),
            TempLabels = [TempLabels handles.SSASA.UnDirFileInfo.NoteLabels{i}];
        end
        TempLabels = unique(TempLabels);
        for i = 1:length(TempLabels),
            handles.SSASA.SyllLabels{i} = TempLabels(i);
        end
        handles.SSASA.SyllCounter = 1;
        set(handles.SyllableListBox, 'String', handles.SSASA.SyllLabels);
        set(handles.SyllableListBox, 'Value', handles.SSASA.SyllCounter);
        
        handles.SSASA.SyllExampleCounter = 1;
        handles.SSASA.SyllData = SSASAGetIndividualSyllData(handles);
        
        % Variables for plotting
        Dir = handles.SSASA.RawDataDirectory;
        FileName = handles.SSASA.SyllData{handles.SSASA.SyllCounter}.FileName{handles.SSASA.SyllExampleCounter};
        FileType = handles.SSASA.FileType;
        SpecAxis = handles.ASSLExampleSpectAxes;
        SpikeTrainAxis = handles.ASSLExampleAmplitudeAxes;
        SpikeTimes = handles.SSASA.SyllData{handles.SSASA.SyllCounter}.SpikeTimes{handles.SSASA.SyllExampleCounter}; 
        TimeRange = handles.SSASA.SyllData{handles.SSASA.SyllCounter}.SyllDur{handles.SSASA.SyllExampleCounter};
        TimeRange(1) = TimeRange(1) - 0.2;
        TimeRange(2) = TimeRange(2) + 0.2;
        FileDur = handles.SSASA.SyllData{handles.SSASA.SyllCounter}.FileDur{handles.SSASA.SyllExampleCounter};
        
        SSAPlotExampleSpectST(Dir, FileName, FileType, SpecAxis, SpikeTrainAxis, SpikeTimes, TimeRange, FileDur, 0.2, 0.2);
        
        axes(handles.ASSLAllSyllWaveformAxes);
        cla(handles.ASSLAllSyllWaveformAxes);
        PlotRaster(handles.SSASA.SyllData{handles.SSASA.SyllCounter}.UnWarpedRaster, 'k');
        axis(handles.SSASA.SyllData{handles.SSASA.SyllCounter}.UnWarpedRasterAxis);
    else
        set(handles.SyllableListBox, 'String', []);
    end
    
end
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SSASyllableAnalysis wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = SSASyllableAnalysis_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in PrevExampleButton.
function PrevExampleButton_Callback(hObject, eventdata, handles)
% hObject    handle to PrevExampleButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.SSASA.SyllExampleCounter = handles.SSASA.SyllExampleCounter - 1;

if (handles.SSASA.SyllExampleCounter < 1)
    handles.SSASA.SyllExampleCounter = 1;
end

if (handles.SSASA.SyllExampleCounter > length(handles.SSASA.SyllData{handles.SSASA.SyllCounter}.FileName))
    handles.SSASA.SyllExampleCounter = length(handles.SSASA.SyllData{handles.SSASA.SyllCounter}.FileName);
end
   
Dir = handles.SSASA.RawDataDirectory;
FileName = handles.SSASA.SyllData{handles.SSASA.SyllCounter}.FileName{handles.SSASA.SyllExampleCounter};
FileType = handles.SSASA.FileType;
SpecAxis = handles.ASSLExampleSpectAxes;
SpikeTrainAxis = handles.ASSLExampleAmplitudeAxes;
SpikeTimes = handles.SSASA.SyllData{handles.SSASA.SyllCounter}.SpikeTimes{handles.SSASA.SyllExampleCounter}; 
TimeRange = handles.SSASA.SyllData{handles.SSASA.SyllCounter}.SyllDur{handles.SSASA.SyllExampleCounter};
TimeRange(1) = TimeRange(1) - 0.2;
TimeRange(2) = TimeRange(2) + 0.2;
FileDur = handles.SSASA.SyllData{handles.SSASA.SyllCounter}.FileDur{handles.SSASA.SyllExampleCounter};

SSAPlotExampleSpectST(Dir, FileName, FileType, SpecAxis, SpikeTrainAxis, SpikeTimes, TimeRange, FileDur, 0.2, 0.2)
guidata(hObject, handles);

% --- Executes on button press in NextExampleButton.
function NextExampleButton_Callback(hObject, eventdata, handles)
% hObject    handle to NextExampleButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.SSASA.SyllExampleCounter = handles.SSASA.SyllExampleCounter + 1;

if (handles.SSASA.SyllExampleCounter < 1)
    handles.SSASA.SyllExampleCounter = 1;
end

if (handles.SSASA.SyllExampleCounter > length(handles.SSASA.SyllData{handles.SSASA.SyllCounter}.FileName))
    handles.SSASA.SyllExampleCounter = length(handles.SSASA.SyllData{handles.SSASA.SyllCounter}.FileName);
end

Dir = handles.SSASA.RawDataDirectory;
FileName = handles.SSASA.SyllData{handles.SSASA.SyllCounter}.FileName{handles.SSASA.SyllExampleCounter};
FileType = handles.SSASA.FileType;
SpecAxis = handles.ASSLExampleSpectAxes;
SpikeTrainAxis = handles.ASSLExampleAmplitudeAxes;
SpikeTimes = handles.SSASA.SyllData{handles.SSASA.SyllCounter}.SpikeTimes{handles.SSASA.SyllExampleCounter}; 
TimeRange = handles.SSASA.SyllData{handles.SSASA.SyllCounter}.SyllDur{handles.SSASA.SyllExampleCounter};
TimeRange(1) = TimeRange(1) - 0.2;
TimeRange(2) = TimeRange(2) + 0.2;
FileDur = handles.SSASA.SyllData{handles.SSASA.SyllCounter}.FileDur{handles.SSASA.SyllExampleCounter};

SSAPlotExampleSpectST(Dir, FileName, FileType, SpecAxis, SpikeTrainAxis, SpikeTimes, TimeRange, FileDur, 0.2, 0.2)
guidata(hObject, handles);

% --- Executes on selection change in SyllableListBox.
function SyllableListBox_Callback(hObject, eventdata, handles)
% hObject    handle to SyllableListBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns SyllableListBox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from SyllableListBox
% Variables for plotting
handles.SSASA.SyllCounter = get(hObject, 'Value');
handles.SSASA.SyllExampleCounter = 1;

Dir = handles.SSASA.RawDataDirectory;
FileName = handles.SSASA.SyllData{handles.SSASA.SyllCounter}.FileName{handles.SSASA.SyllExampleCounter};
FileType = handles.SSASA.FileType;
SpecAxis = handles.ASSLExampleSpectAxes;
SpikeTrainAxis = handles.ASSLExampleAmplitudeAxes;
SpikeTimes = handles.SSASA.SyllData{handles.SSASA.SyllCounter}.SpikeTimes{handles.SSASA.SyllExampleCounter}; 
TimeRange = handles.SSASA.SyllData{handles.SSASA.SyllCounter}.SyllDur{handles.SSASA.SyllExampleCounter};
TimeRange(1) = TimeRange(1) - 0.2;
TimeRange(2) = TimeRange(2) + 0.2;
FileDur = handles.SSASA.SyllData{handles.SSASA.SyllCounter}.FileDur{handles.SSASA.SyllExampleCounter};

SSAPlotExampleSpectST(Dir, FileName, FileType, SpecAxis, SpikeTrainAxis, SpikeTimes, TimeRange, FileDur, 0.2, 0.2);

axes(handles.ASSLAllSyllWaveformAxes);
cla(handles.ASSLAllSyllWaveformAxes);
PlotRaster(handles.SSASA.SyllData{handles.SSASA.SyllCounter}.UnWarpedRaster, 'k');
axis(handles.SSASA.SyllData{handles.SSASA.SyllCounter}.UnWarpedRasterAxis);
        
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function SyllableListBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SyllableListBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
