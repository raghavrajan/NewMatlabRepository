function varargout = DIAPlotRawSpikeData(varargin)
% DIAPLOTRAWSPIKEDATA M-file for DIAPlotRawSpikeData.fig
%      DIAPLOTRAWSPIKEDATA, by itself, creates a new DIAPLOTRAWSPIKEDATA or raises the existing
%      singleton*.
%
%      H = DIAPLOTRAWSPIKEDATA returns the handle to a new DIAPLOTRAWSPIKEDATA or the handle to
%      the existing singleton*.
%
%      DIAPLOTRAWSPIKEDATA('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DIAPLOTRAWSPIKEDATA.M with the given input arguments.
%
%      DIAPLOTRAWSPIKEDATA('Property','Value',...) creates a new DIAPLOTRAWSPIKEDATA or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before DIAPlotRawSpikeData_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to DIAPlotRawSpikeData_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help DIAPlotRawSpikeData

% Last Modified by GUIDE v2.5 25-Aug-2009 12:25:59

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DIAPlotRawSpikeData_OpeningFcn, ...
                   'gui_OutputFcn',  @DIAPlotRawSpikeData_OutputFcn, ...
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


% --- Executes just before DIAPlotRawSpikeData is made visible.
function DIAPlotRawSpikeData_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to DIAPlotRawSpikeData (see VARARGIN)

% Choose default command line output for DIAPlotRawSpikeData
handles.output = hObject;
handles.DIAData = varargin{1};
handles.FileCounter = 1;
handles.ZoomState = 0;

% Update handles structure
guidata(hObject, handles);

DirectoryName = handles.DIAData.DataDirectoryName;
FileName = handles.DIAData.DataFiles(handles.FileCounter,:);
ChanNo = handles.DIAData.SpikeChannelNo;

[Data, Fs] = ReadOKrankData(DirectoryName, FileName, ChanNo);
Data = Data * 100;

Time = (1:1:length(Data))/Fs;
axes(handles.RawSpikeDataPlot);
plot(Time, Data);
hold on;
axis tight;
plot([Time(1) Time(end)], [handles.DIAData.Threshold1 handles.DIAData.Threshold1], 'r');
plot([Time(1) Time(end)], [handles.DIAData.Threshold2 handles.DIAData.Threshold2], 'g');
xlabel('Time (sec)', 'FontWeight', 'bold', 'FontSize', 12);
ylabel('Amplitude (uV)', 'FontWeight', 'bold', 'FontSize', 12);
title(FileName, 'FontWeight', 'bold', 'FontSize', 12);

if (isfield(handles.DIAData, 'Spikes'))
    SpikeTimes = [handles.DIAData.Spikes{handles.FileCounter}];
    for i = 1:length(SpikeTimes),
        WaveformIndices = (round(SpikeTimes(i) * Fs) - 8):1:(round(SpikeTimes(i) * Fs) + 23);
        plot(Time(WaveformIndices), Data(WaveformIndices),'r');
    end
end
    

% UIWAIT makes DIAPlotRawSpikeData wait for user response (see UIRESUME)
% uiwait(handles.DIARawSpikeDataFigure);


% --- Outputs from this function are returned to the command line.
function varargout = DIAPlotRawSpikeData_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in PreviousButton.
function PreviousButton_Callback(hObject, eventdata, handles)
% hObject    handle to PreviousButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (handles.FileCounter == 1)
    return;
end

handles.FileCounter = handles.FileCounter - 1;

% Update handles structure
guidata(hObject, handles);

DirectoryName = handles.DIAData.DataDirectoryName;
FileName = handles.DIAData.DataFiles(handles.FileCounter,:);
ChanNo = handles.DIAData.SpikeChannelNo;

[Data, Fs] = ReadOKrankData(DirectoryName, FileName, ChanNo);
Data = Data * 100;
Time = (1:1:length(Data))/Fs;
axes(handles.RawSpikeDataPlot);
hold off;
plot(Time, Data);
hold on;
axis tight;
plot([Time(1) Time(end)], [handles.DIAData.Threshold1 handles.DIAData.Threshold1], 'r');
plot([Time(1) Time(end)], [handles.DIAData.Threshold2 handles.DIAData.Threshold2], 'g');
xlabel('Time (sec)', 'FontWeight', 'bold', 'FontSize', 12);
ylabel('Amplitude (uV)', 'FontWeight', 'bold', 'FontSize', 12);
title(FileName, 'FontWeight', 'bold', 'FontSize', 12);

if (isfield(handles.DIAData, 'Spikes'))
    SpikeTimes = [handles.DIAData.Spikes{handles.FileCounter}];
    for i = 1:length(SpikeTimes),
        WaveformIndices = (round(SpikeTimes(i) * Fs) - 8):1:(round(SpikeTimes(i) * Fs) + 23);
        plot(Time(WaveformIndices), Data(WaveformIndices),'r');
    end
end

% --- Executes on button press in NextButton.
function NextButton_Callback(hObject, eventdata, handles)
% hObject    handle to NextButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (handles.FileCounter == size(handles.DIAData.DataFiles,1))
    return;
end

handles.FileCounter = handles.FileCounter + 1;

% Update handles structure
guidata(hObject, handles);

DirectoryName = handles.DIAData.DataDirectoryName;
FileName = handles.DIAData.DataFiles(handles.FileCounter,:);
ChanNo = handles.DIAData.SpikeChannelNo;

[Data, Fs] = ReadOKrankData(DirectoryName, FileName, ChanNo);
Data = Data * 100;
Time = (1:1:length(Data))/Fs;
axes(handles.RawSpikeDataPlot);
hold off;
plot(Time, Data);
hold on;
axis tight;
plot([Time(1) Time(end)], [handles.DIAData.Threshold1 handles.DIAData.Threshold1], 'r');
plot([Time(1) Time(end)], [handles.DIAData.Threshold2 handles.DIAData.Threshold2], 'g');
xlabel('Time (sec)', 'FontWeight', 'bold', 'FontSize', 12);
ylabel('Amplitude (uV)', 'FontWeight', 'bold', 'FontSize', 12);
title(FileName, 'FontWeight', 'bold', 'FontSize', 12);

if (isfield(handles.DIAData, 'Spikes'))
    SpikeTimes = [handles.DIAData.Spikes{handles.FileCounter}];
    for i = 1:length(SpikeTimes),
        WaveformIndices = (round(SpikeTimes(i) * Fs) - 8):1:(round(SpikeTimes(i) * Fs) + 23);
        plot(Time(WaveformIndices), Data(WaveformIndices),'r');
    end
end

% --- Executes on button press in CloseButton.
function CloseButton_Callback(hObject, eventdata, handles)
% hObject    handle to CloseButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close(findobj('Tag', 'DIARawSpikeDataFigure'));


% --- Executes on button press in ZoomButton.
function ZoomButton_Callback(hObject, eventdata, handles)
% hObject    handle to ZoomButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (handles.ZoomState == 0)
    handles.ZoomState = 1;
    zoom on;
else
    handles.ZoomState = 0;
    zoom out;
    zoom off;
end
guidata(hObject, handles);