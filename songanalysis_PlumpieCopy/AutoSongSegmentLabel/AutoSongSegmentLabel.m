function varargout = AutoSongSegmentLabel(varargin)
% AUTOSONGSEGMENTLABEL M-file for AutoSongSegmentLabel.fig
%      AUTOSONGSEGMENTLABEL, by itself, creates a new AUTOSONGSEGMENTLABEL or raises the existing
%      singleton*.
%
%      H = AUTOSONGSEGMENTLABEL returns the handle to a new AUTOSONGSEGMENTLABEL or the handle to
%      the existing singleton*.
%
%      AUTOSONGSEGMENTLABEL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in AUTOSONGSEGMENTLABEL.M with the given input arguments.
%
%      AUTOSONGSEGMENTLABEL('Property','Value',...) creates a new AUTOSONGSEGMENTLABEL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before AutoSongSegmentLabel_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to AutoSongSegmentLabel_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help AutoSongSegmentLabel

% Last Modified by GUIDE v2.5 18-Apr-2014 10:07:40

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @AutoSongSegmentLabel_OpeningFcn, ...
                   'gui_OutputFcn',  @AutoSongSegmentLabel_OutputFcn, ...
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


% --- Executes just before AutoSongSegmentLabel is made visible.
function AutoSongSegmentLabel_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to AutoSongSegmentLabel (see VARARGIN)

% Choose default command line output for AutoSongSegmentLabel
handles.output = hObject;

% Variables used by Auto Song Segment and Label
handles.ASSL.FileType = 'okrank';
set(handles.OKrankFileButton, 'Value', 1);
set(handles.ObsFileButton, 'Value', 0);
set(handles.WavFileButton, 'Value', 0);

handles.ASSL.InputFileMode = 'single';
set(handles.SingleFileButton, 'Value', 1);
set(handles.FileListButton, 'Value', 0);

handles.ASSL.FileIndex = 1;

handles.ASSL.SongChanNo = str2double(get(handles.SongChanNoEdit, 'String'));

set(handles.XOffsetSlider, 'Value', 0);
set(handles.XZoomSlider, 'Value', 1);
handles.ASSL.XOffset = str2double(get(handles.XOffsetSlider, 'Value'));
handles.ASSL.XZoom = str2double(get(handles.XZoomSlider, 'Value'));

handles.ASSL.FFTWinSizeSegmenting = str2double(get(handles.FFTWinSizeSegmentingEdit, 'String'));
handles.ASSL.FFTWinOverlapSegmenting = str2double(get(handles.FFTWinOverlapSegmentingEdit, 'String'));
handles.ASSL.MinInt = str2double(get(handles.MinIntEdit, 'String'));
handles.ASSL.MinDur = str2double(get(handles.MinSyllDurEdit, 'String'));

handles.ASSL.FFTWinOverlapTempMatch = str2double(get(handles.FFTWinOverlapTempMatchEdit, 'String'));
handles.ASSL.FFTWinSizeTempMatch = str2double(get(handles.FFTWinSizeTempMatchEdit, 'String'));

handles.ASSL.AutoThreshold = get(handles.AutoThresholdButton, 'Value');
if (handles.ASSL.AutoThreshold == 0)
    set(handles.ThresholdEdit, 'Enable', 'on');
    handles.ASSL.FixedThreshold = get(handles.ThresholdEdit, 'Value');
else
    set(handles.ThresholdEdit, 'Enable', 'off');
end

ProgramDir = which('AutoSongSegmentLabel');
SlashIndex = find((ProgramDir == '/') | (ProgramDir == '\'));
ProgramDir = ProgramDir(1:SlashIndex(end));
Fid = fopen([ProgramDir, 'ASSLFeatureList.txt'], 'r');
TempFeatures = textscan(Fid, '%s', 'DeLimiter', '\n');
fclose(Fid);

handles.ASSL.ToBeUsedFeatures = TempFeatures{1};
set(handles.ToBeUsedFeaturesListBox, 'String', handles.ASSL.ToBeUsedFeatures);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes AutoSongSegmentLabel wait for user response (see UIRESUME)
% uiwait(handles.SpikeSorterMainFig);


% --- Outputs from this function are returned to the command line.
function varargout = AutoSongSegmentLabel_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% ============ Start of Load Files callbacks =============================%

% --- Executes on button press in OKrankFileButton.
function OKrankFileButton_Callback(hObject, eventdata, handles)
% hObject    handle to OKrankFileButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of OKrankFileButton
Val = get(hObject, 'Value');
if (Val == 1)
    handles.ASSL.FileType = 'okrank';
    set(handles.ObsFileButton, 'Value', 0);
    set(handles.WavFileButton, 'Value', 0);
else
    handles.ASSL.FileType = 'obs';
    set(handles.ObsFileButton, 'Value', 1);
    set(handles.WavFileButton, 'Value', 0);
end
guidata(hObject, handles);
% --- Executes on button press in ObsFileButton.
function ObsFileButton_Callback(hObject, eventdata, handles)
% hObject    handle to ObsFileButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ObsFileButton
Val = get(hObject, 'Value');
if (Val == 1)
    handles.ASSL.FileType = 'obs';
    set(handles.OKrankFileButton, 'Value', 0);
    set(handles.WavFileButton, 'Value', 0);
else
    handles.ASSL.FileType = 'wav';
    set(handles.OKrankFileButton, 'Value', 0);
    set(handles.WavFileButton, 'Value', 1);
end
guidata(hObject, handles);

% --- Executes on button press in WavFileButton.
function WavFileButton_Callback(hObject, eventdata, handles)
% hObject    handle to WavFileButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of WavFileButton
Val = get(hObject, 'Value');
if (Val == 1)
    handles.ASSL.FileType = 'wav';
    set(handles.OKrankFileButton, 'Value', 0);
    set(handles.ObsFileButton, 'Value', 0);
else
    handles.ASSL.FileType = 'okrank';
    set(handles.OKrankFileButton, 'Value', 1);
    set(handles.ObsFileButton, 'Value', 0);
end
guidata(hObject, handles);

% --- Executes on button press in SingleFileButton.
function SingleFileButton_Callback(hObject, eventdata, handles)
% hObject    handle to SingleFileButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of SingleFileButton
Val = get(hObject, 'Value');
if (Val == 1)
    handles.ASSL.InputFileMode = 'single';
    set(handles.FileListButton, 'Value', 0);
else
    handles.ASSL.InputFileMode = 'multiple';
    set(handles.FileListButton, 'Value', 1);
end
guidata(hObject, handles);

% --- Executes on button press in FileListButton.
function FileListButton_Callback(hObject, eventdata, handles)
% hObject    handle to FileListButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Val = get(hObject, 'Value');
if (Val == 1)
    handles.ASSL.InputFileMode = 'multiple';
    set(handles.SingleFileButton, 'Value', 0);
else
    handles.ASSL.InputFileMode = 'single';
    set(handles.SingleFileButton, 'Value', 1);
end
guidata(hObject, handles);
% Hint: get(hObject,'Value') returns toggle state of FileListButton


% --- Executes on button press in LoadFilesButton.
function [RawData, Fs, Time] = LoadFilesButton_Callback(hObject, eventdata, handles, varargin)
% hObject    handle to LoadFilesButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isempty(strfind(handles.ASSL.InputFileMode, 'single'))
    if (nargin < 4)
        [handles.ASSL.FileName{1}, handles.ASSL.DirName] = uigetfile('*', 'Choose input file');
        if (isfield(handles.ASSL, 'Threshold'))
            handles = rmfield(handles.ASSL, 'Threshold');
        end
        if (isfield(handles.ASSL, 'SyllOnsets'))
            handles = rmfield(handles.ASSL, 'SyllOnsets');
            handles = rmfield(handles.ASSL, 'SyllOffsets');
        end
        if (isfield(handles.ASSL, 'SyllLabels'))
            handles = rmfield(handles.ASSL, 'SyllLabels');
        end
    end
    if (ispc)
        if (handles.ASSL.DirName(end) ~= '\')
            handles.ASSL.DirName(end+1) = '\';
        else
            handles.ASSL.DirName(end+1) = '/';
        end
    end
    handles.ASSL.FileIndex = 1;
else
    if (nargin < 4)
        [handles.ASSL.FileListName, handles.ASSL.DirName] = uigetfile('*', 'Choose input file list');
    end
    
    if (ispc)
        if (handles.ASSL.DirName(end) ~= '\')
            handles.ASSL.DirName(end+1) = '\';
        else
            handles.ASSL.DirName(end+1) = '/';
        end
    end
    
    Fid = fopen([handles.ASSL.DirName, handles.ASSL.FileListName], 'r');
    TempFiles = textscan(Fid, '%s', 'DeLimiter', '\n');
    handles.ASSL.FileName = TempFiles{1};
    fclose(Fid);

    for i = 1:length(handles.ASSL.FileName),
        SlashIndex = find((handles.ASSL.FileName{i} == '/') | (handles.ASSL.FileName{i} == '\'));
        if (~isempty(SlashIndex))
            handles.ASSL.FileName{i} = handles.ASSL.FileName{i}(SlashIndex(end)+1:end);
        end
    end
end

[RawData, Fs] = ASSLGetRawData(handles.ASSL.DirName, handles.ASSL.FileName{handles.ASSL.FileIndex}, handles.ASSL.FileType, handles.ASSL.SongChanNo);

Time = (1:1:length(RawData))/Fs;

[LogAmplitude] = ASSLCalculateLogAmplitude(RawData, Fs, Time, handles.ASSL.FFTWinSizeSegmenting, handles.ASSL.FFTWinOverlapSegmenting);

set(handles.SongFileNameText, 'String', ['Song File Name : ', handles.ASSL.FileName{handles.ASSL.FileIndex}]);

if (isfield(handles.ASSL, 'SyllOnsets'))
    if (isfield(handles.ASSL, 'SyllLabels'))
        [handles.ASSL.SpecAxisLimits, handles.ASSL.AmpAxisLimits] = ASSLPlotData(handles.ASSL.DirName, handles.ASSL.FileName{handles.ASSL.FileIndex}, handles.ASSL.FileType, Time, LogAmplitude, handles.ASSLSpecAxis, handles.ASSLAmpAxis, handles.ASSL.Threshold{handles.ASSL.FileIndex}, handles.ASSL.SyllOnsets{handles.ASSL.FileIndex}, handles.ASSL.SyllOffsets{handles.ASSL.FileIndex}, handles.ASSL.SyllLabels{handles.ASSL.FileIndex});
    else
        [handles.ASSL.SpecAxisLimits, handles.ASSL.AmpAxisLimits] = ASSLPlotData(handles.ASSL.DirName, handles.ASSL.FileName{handles.ASSL.FileIndex}, handles.ASSL.FileType, Time, LogAmplitude, handles.ASSLSpecAxis, handles.ASSLAmpAxis, handles.ASSL.Threshold{handles.ASSL.FileIndex}, handles.ASSL.SyllOnsets{handles.ASSL.FileIndex}, handles.ASSL.SyllOffsets{handles.ASSL.FileIndex});
    end
else
    [handles.ASSL.SpecAxisLimits, handles.ASSL.AmpAxisLimits] = ASSLPlotData(handles.ASSL.DirName, handles.ASSL.FileName{handles.ASSL.FileIndex}, handles.ASSL.FileType, Time, LogAmplitude, handles.ASSLSpecAxis, handles.ASSLAmpAxis);
end

set(handles.XOffsetSlider, 'Max', Time(end));
set(handles.XOffsetSlider, 'Min', Time(1));
set(handles.XOffsetSlider, 'Value', Time(1));

set(handles.XZoomSlider, 'Max', Time(end));
set(handles.XZoomSlider, 'Min', 0.01);
set(handles.XZoomSlider, 'Value', get(handles.XZoomSlider, 'Max'));

handles.ASSL.XOffset = get(handles.XOffsetSlider, 'Value');
handles.ASSL.XZoom = get(handles.XZoomSlider, 'Value');


FileSep = filesep;

if (handles.ASSL.DirName(end) ~= FileSep)
    handles.ASSL.NoteFileDirName = [handles.ASSL.DirName, FileSep, 'ASSLNoteFiles', FileSep];
else
    handles.ASSL.NoteFileDirName = [handles.ASSL.DirName, 'ASSLNoteFiles', FileSep];
end

if (~exist(handles.ASSL.NoteFileDirName, 'dir'))
    mkdir(handles.ASSL.NoteFileDirName);
end

guidata(hObject, handles);
        

function SongChanNoEdit_Callback(hObject, eventdata, handles)
% hObject    handle to SongChanNoEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SongChanNoEdit as text
%        str2double(get(hObject,'String')) returns contents of SongChanNoEdit as a double
if (isnan(str2double(get(hObject, 'String'))))
    handles.ASSL.SongChanNo = get(hObject, 'String');
else
    handles.ASSL.SongChanNo = str2double(get(hObject, 'String'));
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function SongChanNoEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SongChanNoEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% ============ End of Load Files callbacks ===============================%

% ============ Start of Syllable Detection callbacks =====================%

% --- Executes on button press in SegmentSongsButton.
function SegmentSongsButton_Callback(hObject, eventdata, handles)
% hObject    handle to SegmentSongsButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


for i = 1:length(handles.ASSL.FileName),
    if (exist([handles.ASSL.NoteFileDirName, handles.ASSL.FileName{i}, '.not.mat'], 'file'))
        Temp = load([handles.ASSL.NoteFileDirName, handles.ASSL.FileName{i}, '.not.mat']);
        handles.ASSL.SyllOnsets{i} = Temp.onsets;
        handles.ASSL.SyllOffsets{i} = Temp.offsets;
        handles.ASSL.SyllLabels{i} = Temp.labels;
        if (isfield(Temp, 'threshold'))
            handles.ASSL.Threshold{i} = Temp.threshold;
        else
            handles.ASSL.Threshold{i} = -30;
        end
        disp(['Loaded syllable data from ', handles.ASSL.FileName{i}, '.not.mat']);
        if (~isfield(Temp, 'FileLength'))
            [RawData, Fs] = ASSLGetRawData(handles.ASSL.DirName, handles.ASSL.FileName{i}, handles.ASSL.FileType, handles.ASSL.SongChanNo);
            handles.ASSL.FileDur{i} = length(RawData)/Fs;
        end
        continue;
    end
    [RawData, Fs] = ASSLGetRawData(handles.ASSL.DirName, handles.ASSL.FileName{i}, handles.ASSL.FileType, handles.ASSL.SongChanNo);
    handles.ASSL.FileDur{i} = length(RawData)/Fs;
    
    Time = (1:1:length(RawData))/Fs;
    [LogAmplitude] = ASSLCalculateLogAmplitude(RawData, Fs, Time, handles.ASSL.FFTWinSizeSegmenting, handles.ASSL.FFTWinOverlapSegmenting);
    
    if (handles.ASSL.AutoThreshold == 0)
        handles.ASSL.Threshold{i} = handles.ASSL.FixedThreshold;
    else
        %Obj = gmdistribution.fit(LogAmplitude', 2);
        handles.ASSL.Threshold{i} = ASSLCalculateFisherThreshold(LogAmplitude);
        %handles.ASSL.Threshold{i} = mean(Obj.mu);
    end
    
    [handles.ASSL.SyllOnsets{i}, handles.ASSL.SyllOffsets{i}] = ASSLSegmentData(LogAmplitude, Fs, handles.ASSL.MinInt, handles.ASSL.MinDur, handles.ASSL.Threshold{i});
    for j = 1:length(handles.ASSL.SyllOnsets{i}),
        handles.ASSL.SyllLabels{i}(j) = '0';
    end
    disp([' Detected ', num2str(length(handles.ASSL.SyllOnsets{i})), ' syllables for ', handles.ASSL.FileName{i}]);
end

set(handles.SongFileNameText, 'String', ['Song File Name : ', handles.ASSL.FileName{handles.ASSL.FileIndex}]);
[RawData, Fs] = ASSLGetRawData(handles.ASSL.DirName, handles.ASSL.FileName{handles.ASSL.FileIndex}, handles.ASSL.FileType, handles.ASSL.SongChanNo);
    
Time = (1:1:length(RawData))/Fs;
[LogAmplitude] = ASSLCalculateLogAmplitude(RawData, Fs, Time, handles.ASSL.FFTWinSizeSegmenting, handles.ASSL.FFTWinOverlapSegmenting);

[handles.ASSL.SpecAxisLimits, handles.ASSL.AmpAxisLimits] = ASSLPlotData(handles.ASSL.DirName, handles.ASSL.FileName{handles.ASSL.FileIndex}, handles.ASSL.FileType, Time, LogAmplitude, handles.ASSLSpecAxis, handles.ASSLAmpAxis, handles.ASSL.Threshold{handles.ASSL.FileIndex}, handles.ASSL.SyllOnsets{handles.ASSL.FileIndex}, handles.ASSL.SyllOffsets{handles.ASSL.FileIndex});

set(handles.XOffsetSlider, 'Max', Time(end));
set(handles.XZoomSlider, 'Max', Time(end));
set(handles.XZoomSlider, 'Min', 0.01);
set(handles.XZoomSlider, 'Value', get(handles.XZoomSlider, 'Max'));
    
guidata(hObject, handles);

function FFTWinSizeSegmentingEdit_Callback(hObject, eventdata, handles)
% hObject    handle to FFTWinSizeSegmentingEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FFTWinSizeSegmentingEdit as text
%        str2double(get(hObject,'String')) returns contents of FFTWinSizeSegmentingEdit as a double
handles.ASSL.FFTWinSizeSegmenting = str2double(get(hObject, 'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function FFTWinSizeSegmentingEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FFTWinSizeSegmentingEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function MinIntEdit_Callback(hObject, eventdata, handles)
% hObject    handle to MinIntEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MinIntEdit as text
%        str2double(get(hObject,'String')) returns contents of MinIntEdit as a double

handles.ASSL.MinInt = str2double(get(hObject, 'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function MinIntEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MinIntEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function FFTWinOverlapSegmentingEdit_Callback(hObject, eventdata, handles)
% hObject    handle to FFTWinOverlapSegmentingEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FFTWinOverlapSegmentingEdit as text
%        str2double(get(hObject,'String')) returns contents of FFTWinOverlapSegmentingEdit as a double

handles.ASSL.FFTWinOverlapSegmenting = str2double(get(hObject, 'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function FFTWinOverlapSegmentingEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FFTWinOverlapSegmentingEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in CalculateFeaturesButton.
function CalculateFeaturesButton_Callback(hObject, eventdata, handles)
% hObject    handle to CalculateFeaturesButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.ASSL.FeatValues = [];
handles.ASSL.RawFeatValues = [];
handles.ASSL.SyllIndices = [];
handles.ASSL.SyllIndexLabels = [];

Filesep = filesep;

for i = 1:length(handles.ASSL.ToBeUsedFeatures),
    if (handles.ASSL.NoteFileDirName(end) == Filesep)
        if ((isfield(handles.ASSL, 'FileListName')) && (~isempty(handles.ASSL.FileListName)))
            RootOutputFileName = [handles.ASSL.NoteFileDirName, handles.ASSL.FileListName];
        else
            RootOutputFileName = [handles.ASSL.NoteFileDirName, handles.ASSL.FileName{1}];
        end
    else
        if ((isfield(handles.ASSL, 'FileListName')) && (~isempty(handles.ASSL.FileListName)))
            RootOutputFileName = [handles.ASSL.NoteFileDirName, Filesep, handles.ASSL.FileListName];
        else
            RootOutputFileName = [handles.ASSL.NoteFileDirName, Filesep, handles.ASSL.FileName{1}];
        end
    end
    
    FeatValuesFileName = [RootOutputFileName, '.FeatValues.mat'];
    RawFeatValuesFileName = [RootOutputFileName, '.RawFeatValues.mat'];
    SyllIndicesFileName = [RootOutputFileName, '.SyllIndices.mat'];
    SyllIndexLabelsFileName = [RootOutputFileName, '.SyllIndexLabels.mat'];
    RawOutputFileName = [RootOutputFileName, '.', handles.ASSL.ToBeUsedFeatures{i}, '.raw.mat'];
    OutputFileName = [RootOutputFileName, '.', handles.ASSL.ToBeUsedFeatures{i}, '.mat'];
    
    if (exist(OutputFileName, 'file'))
        load(OutputFileName);
        eval(['handles.ASSL.', handles.ASSL.ToBeUsedFeatures{i}, ' = Temp;']);
        if (((isempty(strfind(handles.ASSL.ToBeUsedFeatures{i}, 'Duration')))) && ((isempty(strfind(handles.ASSL.ToBeUsedFeatures{i}, 'EntropyVariance')))))
            if (exist(RawOutputFileName, 'file'))
                load(RawOutputFileName);
                eval(['handles.ASSL.Raw.', handles.ASSL.ToBeUsedFeatures{i}, ' = Temp;']);
                Flag(i) = 1;
            else
                Flag(i) = 0;
            end
        else
            if (exist(FeatValuesFileName, 'file'))
                load(FeatValuesFileName);
                eval(['handles.ASSL.FeatValues = Temp;']);
                Flag(i) = 1;
            else
                Flag(i) = 0;
            end
            if (exist(RawFeatValuesFileName, 'file'))
                load(RawFeatValuesFileName);
                eval(['handles.ASSL.RawFeatValues = Temp;']);
            end
            if (exist(SyllIndicesFileName, 'file'))
                load(SyllIndicesFileName);
                eval(['handles.ASSL.SyllIndices = Temp;']);
            end
            if (exist(SyllIndexLabelsFileName, 'file'))
                load(SyllIndexLabelsFileName);
                eval(['handles.ASSL.SyllIndexLabels = Temp;']);
            end
        end
    else
        Flag(i) = 0;
    end
end

if (sum(Flag) == length(handles.ASSL.ToBeUsedFeatures))
    disp('Loaded calculated feature values from existing files');
    set(handles.SaveDataButton, 'enable', 'on');
    guidata(hObject, handles);
    return;
end

handles.ASSL.FeatValues = [];
handles.ASSL.RawFeatValues = [];
handles.ASSL.SyllIndices = [];
handles.ASSL.SyllIndexLabels = [];

SyllNo = 0;
fprintf('\n');
for i = 1:length(handles.ASSL.FileName),
    if (mod(i,5) == 0),
        fprintf('>');
    end
    TempFeats = [];
    TempRawFeats = [];
    [RawData, Fs] = ASSLGetRawData(handles.ASSL.DirName, handles.ASSL.FileName{i}, handles.ASSL.FileType, handles.ASSL.SongChanNo);
    
    Time = (1:1:length(RawData))/Fs;
    
    [Feats, RawFeats] = ASSLCalculateSAPFeatsWithOnsets(RawData, Time, Fs, handles.ASSL.SyllOnsets{i}/1000, handles.ASSL.SyllOffsets{i}/1000);
    
    for j = 1:length(handles.ASSL.ToBeUsedFeatures),
        FeatFields = fieldnames(Feats);
        FeatureIndex = find(strcmp(FeatFields, handles.ASSL.ToBeUsedFeatures{j}));
        eval(['handles.ASSL.', FeatFields{FeatureIndex}, '{', num2str(i), '} = Feats.', FeatFields{FeatureIndex}, ';']);
        TempFeats = [TempFeats eval(['Feats.', FeatFields{FeatureIndex}])'];
        if (isfield(RawFeats, handles.ASSL.ToBeUsedFeatures{j}))
            eval(['handles.ASSL.Raw.', FeatFields{FeatureIndex}, '{', num2str(i), '} = RawFeats.', FeatFields{FeatureIndex}, ';']);
            TempRawFeats = [TempRawFeats eval(['RawFeats.', FeatFields{FeatureIndex}])'];
        end
    end
    handles.ASSL.FeatValues = [handles.ASSL.FeatValues; TempFeats];
    handles.ASSL.RawFeatValues = [handles.ASSL.RawFeatValues; TempRawFeats];
    handles.ASSL.SyllIndices = [handles.ASSL.SyllIndices; ones(length(handles.ASSL.SyllOnsets{i}),1)*i [1:1:length(handles.ASSL.SyllOnsets{i})]' ([1:1:length(handles.ASSL.SyllOnsets{i})]' + SyllNo)];
    if (isfield(handles.ASSL, 'SyllLabels'))
        handles.ASSL.SyllIndexLabels = [handles.ASSL.SyllIndexLabels; handles.ASSL.SyllLabels{i}'];
    end
    SyllNo = SyllNo + length(handles.ASSL.SyllOnsets{i});
end
fprintf('\n');
set(handles.SaveDataButton, 'enable', 'on');

guidata(hObject, handles);
disp('Finished calculating features');

% --- Executes on button press in SaveFeatValsButton.
function SaveFeatValsButton_Callback(hObject, eventdata, handles)
% hObject    handle to SaveFeatValsButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Filesep = filesep;

for i = 1:length(handles.ASSL.ToBeUsedFeatures),
    if (handles.ASSL.NoteFileDirName(end) == Filesep)
        if ((isfield(handles.ASSL, 'FileListName')) && (~isempty(handles.ASSL.FileListName)))
            RootOutputFileName = [handles.ASSL.NoteFileDirName, handles.ASSL.FileListName];
        else
            RootOutputFileName = [handles.ASSL.NoteFileDirName, handles.ASSL.FileName{1}];
        end
    else
        if ((isfield(handles.ASSL, 'FileListName')) && (~isempty(handles.ASSL.FileListName)))
            RootOutputFileName = [handles.ASSL.NoteFileDirName, Filesep, handles.ASSL.FileListName];
        else
            RootOutputFileName = [handles.ASSL.NoteFileDirName, Filesep, handles.ASSL.FileName{1}];
        end
    end
    
    if (i == 1)
        FeatValuesFileName = [RootOutputFileName, '.FeatValues.mat'];
        Temp = eval(['handles.ASSL.FeatValues']);
        save(FeatValuesFileName, 'Temp');

        RawFeatValuesFileName = [RootOutputFileName, '.RawFeatValues.mat'];
        Temp = eval(['handles.ASSL.RawFeatValues']);
        save(RawFeatValuesFileName, 'Temp');

        SyllIndicesFileName = [RootOutputFileName, '.SyllIndices.mat'];
        Temp = eval(['handles.ASSL.SyllIndices']);
        save(SyllIndicesFileName, 'Temp');

        SyllIndexLabelsFileName = [RootOutputFileName, '.SyllIndexLabels.mat'];
        Temp = eval(['handles.ASSL.SyllIndexLabels']);
        save(SyllIndexLabelsFileName, 'Temp');
    end
    
    RawOutputFileName = [RootOutputFileName,'.', handles.ASSL.ToBeUsedFeatures{i}, '.raw.mat'];
    if (((isempty(strfind(handles.ASSL.ToBeUsedFeatures{i}, 'Duration')))) && ((isempty(strfind(handles.ASSL.ToBeUsedFeatures{i}, 'EntropyVariance')))))
        Temp = eval(['handles.ASSL.Raw.', handles.ASSL.ToBeUsedFeatures{i}]);
        save(RawOutputFileName, 'Temp');
    end
    
    OutputFileName = [RootOutputFileName, '.', handles.ASSL.ToBeUsedFeatures{i}, '.mat'];
    Temp = eval(['handles.ASSL.', handles.ASSL.ToBeUsedFeatures{i}]);
    save(OutputFileName, 'Temp');
end
disp('Finished saving feature value files');
guidata(hObject, handles);

% --- Executes on button press in PlotFeaturesButton.
function PlotFeaturesButton_Callback(hObject, eventdata, handles)
% hObject    handle to PlotFeaturesButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

TempFields = handles.ASSL.ToBeUsedFeatures;
ASSLPlotFeatures(handles.ASSL.FeatValues, TempFields, handles.ASSL);


% --- Executes on button press in LoadTemplatesButton.
function LoadTemplatesButton_Callback(hObject, eventdata, handles)
% hObject    handle to LoadTemplatesButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
PresentDir = pwd;

[handles.ASSL.TemplateFileName, handles.ASSL.TemplateDir] = uigetfile('*', 'Choose input file');

cd(handles.ASSL.TemplateDir);

handles.ASSL.Templates{1} = load(handles.ASSL.TemplateFileName);

cd(PresentDir);
guidata(hObject, handles);
disp(['Finished loading ', num2str(length(handles.ASSL.TemplateFileName)), ' templates']);

% --- Executes on button press in DoTemplateMatchingButton.
function DoTemplateMatchingButton_Callback(hObject, eventdata, handles)
% hObject    handle to DoTemplateMatchingButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.ASSL.StretchValues = [-3:3:3];
SyllableIndex = 0;
fprintf('\n');
for i = 1:length(handles.ASSL.FileName),
    if (isempty(handles.ASSL.SyllOnsets{i}))
        handles.ASSL.SyllLabels{i} = [];
        continue;
    end
    handles.ASSL.SyllableTemplateMatches{i} = [];

    fprintf('%s:', handles.ASSL.FileName{i}); 
    [RawData, Fs] = ASSLGetRawData(handles.ASSL.DirName, handles.ASSL.FileName{i}, handles.ASSL.FileType, handles.ASSL.SongChanNo);
    
    Time = (1:1:length(RawData))/Fs;
    
    WinSize = round(handles.ASSL.FFTWinSizeTempMatch/1000 * Fs);
    WinOverlap = round(handles.ASSL.FFTWinOverlapTempMatch * WinSize);

   % [S1, F, T, P] = spectrogram(RawData, hamming(WinSize), WinOverlap, WinSize, Fs);
    [P, F, S1, T] = CalculateMultiTaperSpectrogram(RawData, Fs, 8, 8/2, 1.5);
    
    Freq1 = find((F >= 860) & (F <= 8600));
    S = log10(abs(S1(Freq1,:)));

    %Syllables = find(handles.ASSL.SyllIndices(:,1) == i);
    
    clear TempMatch SyllableMatch MatchLen;
    for k = 1:length(handles.ASSL.Templates{1}.SyllableTemplates),
        
        TempSyllTemplate = handles.ASSL.Templates{1}.SyllableTemplates{k}{min(length(handles.ASSL.Templates{1}.SyllableTemplates{k}),3)};
        
        fprintf('>');
        [TempMatch{k}] = ASSLTemplateMatch(S, TempSyllTemplate.MotifTemplate,  handles.ASSL.StretchValues);
        MatchLen(k) = length(TempMatch{k});
    end
    for k = 1:length(TempMatch),
        SyllableMatch(k,:) = TempMatch{k}(1:min(MatchLen));
    end
    
    SyllableMatchFs = 1/(T(2) - T(1));
    
    for j = 1:length(handles.ASSL.SyllOnsets{i}),
        SyllableOnsetTimeIndex = round(handles.ASSL.SyllOnsets{i}(j)/1000 * SyllableMatchFs);
        SyllableOnsetWindow = [(SyllableOnsetTimeIndex - round(0.015 * SyllableMatchFs)) (SyllableOnsetTimeIndex + round(0.015 * SyllableMatchFs))];
        if (SyllableOnsetWindow(1) <= 0)
            SyllableOnsetWindow(1) = 1;
        end
        
        if ((SyllableOnsetWindow(1) > size(SyllableMatch,2)) || (SyllableOnsetWindow(2) > size(SyllableMatch,2)))
            handles.ASSL.SyllLabels{i}(j) = '0';
            handles.ASSL.SyllableTemplateMatches{i}(j,:) = zeros(1, length(handles.ASSL.Templates{1}.SyllableTemplates));
        else
            [MaxVal, MaxIndex] = max(SyllableMatch(:, SyllableOnsetWindow(1):SyllableOnsetWindow(2)), [], 2);
            handles.ASSL.SyllableTemplateMatches{i}(j,:) = MaxVal;
            [MaxVal, MaxIndex2] = max(MaxVal);
            if (MaxVal > 1)
                handles.ASSL.SyllLabels{i}(j) = handles.ASSL.Templates{1}.SyllableTemplates{MaxIndex2}{1}.MotifTemplate(1).Label;
                TemplateLen = size(handles.ASSL.Templates{1}.SyllableTemplates{MaxIndex2}{1}.MotifTemplate(1).MotifTemplate, 2)/SyllableMatchFs;
                [MinVal, MinIndex] = min(abs(handles.ASSL.SyllOffsets{i}(j:end)/1000 - (handles.ASSL.SyllOnsets{i}(j)/1000 + TemplateLen)));
                handles.ASSL.SyllOffsets{i}(j) = handles.ASSL.SyllOffsets{i}(MinIndex + j - 1);
            else
                handles.ASSL.SyllLabels{i}(j) = '0';
            end
        end
    end
    SyllableCheckIndex = 1;
    while (SyllableCheckIndex < (length(handles.ASSL.SyllOnsets{i})-1))
        Flag = 1;
        for j = 1:length(handles.ASSL.SyllOnsets{i})-1,
            if (Flag == 0)
                break;
            end
            SyllableCheckIndex = j;
            for k = j+1:length(handles.ASSL.SyllOnsets{i}),
                if ((handles.ASSL.SyllOnsets{i}(k) >= handles.ASSL.SyllOnsets{i}(j)) && (handles.ASSL.SyllOnsets{i}(k) <= handles.ASSL.SyllOffsets{i}(j)))
                    handles.ASSL.SyllOnsets{i}(k) = [];
                    handles.ASSL.SyllOffsets{i}(k) = [];
                    handles.ASSL.SyllLabels{i}(k) = [];
                    handles.ASSL.SyllableTemplateMatches{i}(k,:) = [];
                    Flag = 0;
                    break;
                end
            end
        end
    end
    SyllableIndex = SyllableIndex + length(handles.ASSL.SyllOnsets{i});
    fprintf('\n');
end    

set(handles.SongFileNameText, 'String', ['Song File Name : ', handles.ASSL.FileName{handles.ASSL.FileIndex}]);
[RawData, Fs] = ASSLGetRawData(handles.ASSL.DirName, handles.ASSL.FileName{handles.ASSL.FileIndex}, handles.ASSL.FileType, handles.ASSL.SongChanNo);

Time = (1:1:length(RawData))/Fs;
[LogAmplitude] = ASSLCalculateLogAmplitude(RawData, Fs, Time, handles.ASSL.FFTWinSizeSegmenting, handles.ASSL.FFTWinOverlapSegmenting);

[handles.ASSL.SpecAxisLimits, handles.ASSL.AmpAxisLimits] = ASSLPlotData(handles.ASSL.DirName, handles.ASSL.FileName{handles.ASSL.FileIndex}, handles.ASSL.FileType, Time, LogAmplitude, handles.ASSLSpecAxis, handles.ASSLAmpAxis, handles.ASSL.Threshold{handles.ASSL.FileIndex}, handles.ASSL.SyllOnsets{handles.ASSL.FileIndex}, handles.ASSL.SyllOffsets{handles.ASSL.FileIndex}, handles.ASSL.SyllLabels{handles.ASSL.FileIndex});

set(handles.XOffsetSlider, 'Max', Time(end));
set(handles.XZoomSlider, 'Max', Time(end));
set(handles.XZoomSlider, 'Min', 0.01);
set(handles.XZoomSlider, 'Value', get(handles.XZoomSlider, 'Max'));

guidata(hObject, handles);

% --- Executes on button press in ShowResultsButton.
function ShowResultsButton_Callback(hObject, eventdata, handles)
% hObject    handle to ShowResultsButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

TempFields = handles.ASSL.ToBeUsedFeatures;
ASSLTemplateMatchPlotResults(handles.ASSL.FeatValues, TempFields, handles.ASSL);

function MinSyllDurEdit_Callback(hObject, eventdata, handles)
% hObject    handle to MinSyllDurEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MinSyllDurEdit as text
%        str2double(get(hObject,'String')) returns contents of MinSyllDurEdit as a double

handles.ASSL.MinDur = str2double(get(hObject, 'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function MinSyllDurEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MinSyllDurEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function FFTWinOverlapTempMatchEdit_Callback(hObject, eventdata, handles)
% hObject    handle to FFTWinOverlapTempMatchEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FFTWinOverlapTempMatchEdit as text
%        str2double(get(hObject,'String')) returns contents of FFTWinOverlapTempMatchEdit as a double

handles.ASSL.FFTWinOverlapTempMatch = str2double(get(hObject, 'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function FFTWinOverlapTempMatchEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FFTWinOverlapTempMatchEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in WriteNoteFilesButton.
function WriteNoteFilesButton_Callback(hObject, eventdata, handles)
% hObject    handle to WriteNoteFilesButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

for i = 1:length(handles.ASSL.FileName),
    onsets = handles.ASSL.SyllOnsets{i};
    offsets = handles.ASSL.SyllOffsets{i};
    labels = handles.ASSL.SyllLabels{i};
    threshold = handles.ASSL.Threshold{i};
    sm_win = handles.ASSL.FFTWinSizeSegmenting;
    min_int = handles.ASSL.MinInt;
    min_dur = handles.ASSL.MinDur;
    save([handles.ASSL.NoteFileDirName, handles.ASSL.FileName{i}, '.not.mat'], 'onsets', 'offsets', 'labels', 'threshold', 'sm_win', 'min_dur', 'min_int');
    clear onsets offsets labels threshold sm_win min_int min_dur;
end

msgbox('Finished writing syllable boundaries to note files');


% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in WriteAdjustedNoteFilesButton.
function WriteAdjustedNoteFilesButton_Callback(hObject, eventdata, handles)
% hObject    handle to WriteAdjustedNoteFilesButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (~ispc)
    handles.ASB.AdjustDirName = [handles.ASSL.DirName, '/AdjustedNoteFiles/'];
else
    handles.ASB.AdjustDirName = [handles.ASSL.DirName, '\AdjustedNoteFiles\'];
end

if (~exist(handles.ASB.AdjustDirName, 'dir'))
    mkdir(handles.ASB.AdjustDirName);
end

for i = 1:length(handles.ASSL.FileName),
    if (isempty(handles.ASSL.SyllOnsets{i}))
        onsets = [];
        offsets = [];
        labels = [];
    else
        onsets = handles.ASSL.AdjSyllOnsets{i};
        offsets = handles.ASSL.AdjSyllOffsets{i};
        labels = handles.ASSL.SyllLabels{i};
    end
    threshold = handles.ASSL.Threshold{i};
    sm_win = handles.ASSL.FFTWinSizeSegmenting;
    min_int = handles.ASSL.MinInt;
    min_dur = handles.ASSL.MinDur;
    save([handles.ASB.AdjustDirName, handles.ASSL.FileName{i}, '.not.mat'], 'onsets', 'offsets', 'labels', 'threshold', 'sm_win', 'min_dur', 'min_int');
    clear onsets offsets labels threshold sm_win min_int min_dur;
end

msgbox('Finished writing adjusted syllable boundaries to note files');

% --- Executes on button press in LoadDataButton.
function LoadDataButton_Callback(hObject, eventdata, handles)
% hObject    handle to LoadDataButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[SavedDataFile, PathName] = uigetfile('*.mat', 'Choose the file that has the saved data');

load([PathName, SavedDataFile]);
guidata(hObject, handles);
disp(['Loaded data from ', SavedDataFile]);

% --- Executes on button press in WriteLogFileButton.
function WriteLogFileButton_Callback(hObject, eventdata, handles)
% hObject    handle to WriteLogFileButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Time = clock;
OutputFileName = [handles.SS.FileName{1}, '_SSOutput_', num2str(Time(1)), num2str(Time(2)), num2str(Time(3)), num2str(Time(4)), num2str(Time(5)), num2str(round(Time(6))), '.log'];
Fid = fopen(OutputFileName, 'w');


fclose(Fid);

disp(['Wrote log file ', OutputFileName]);


% --- Executes on slider movement.
function XZoomSlider_Callback(hObject, eventdata, handles)
% hObject    handle to XZoomSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.ASSL.XZoom = get(hObject, 'Value');

if (handles.ASSL.XZoom < get(hObject, 'Min'))
    handles.ASSL.XZoom = get(hObject, 'Min');
end

if (handles.ASSL.XZoom > get(hObject, 'Max'))
    handles.ASSL.XZoom = get(hObject, 'Max');
end

axes(handles.ASSLSpecAxis);
handles.ASSL.SpecAxisLimits = [handles.ASSL.XOffset (handles.ASSL.XOffset + handles.ASSL.XZoom) handles.ASSL.SpecAxisLimits(3:4)];
axis(handles.ASSL.SpecAxisLimits);

axes(handles.ASSLAmpAxis);
handles.ASSL.AmpAxisLimits = [handles.ASSL.XOffset (handles.ASSL.XOffset + handles.ASSL.XZoom) handles.ASSL.AmpAxisLimits(3:4)];
axis(handles.ASSL.AmpAxisLimits);

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function XZoomSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to XZoomSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in PrevFileButton.
function PrevFileButton_Callback(hObject, eventdata, handles)
% hObject    handle to PrevFileButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.ASSL.FileIndex = handles.ASSL.FileIndex - 1;
if (handles.ASSL.FileIndex < 1)
    handles.ASSL.FileIndex = 1;
end

if (handles.ASSL.FileIndex > length(handles.ASSL.FileName))
    handles.ASSL.FileIndex = length(handles.ASSL.FileName);
end
guidata(hObject, handles);

LoadFilesButton_Callback(hObject, eventdata, handles, 'PrevFile');

% --- Executes on button press in NextFileButton.
function NextFileButton_Callback(hObject, eventdata, handles)
% hObject    handle to NextFileButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.ASSL.FileIndex = handles.ASSL.FileIndex + 1;

if (handles.ASSL.FileIndex < 1)
    handles.ASSL.FileIndex = 1;
end

if (handles.ASSL.FileIndex > length(handles.ASSL.FileName))
    handles.ASSL.FileIndex = length(handles.ASSL.FileName);
end

guidata(hObject, handles);
LoadFilesButton_Callback(hObject, eventdata, handles, 'NextFile');

% --- Executes on slider movement.
function XOffsetSlider_Callback(hObject, eventdata, handles)
% hObject    handle to XOffsetSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.ASSL.XOffset = get(hObject, 'Value');

if (handles.ASSL.XOffset < get(hObject, 'Min'))
    handles.ASSL.XOffset = get(hObject, 'Min');
end

if (handles.ASSL.XOffset > get(hObject, 'Max'))
    handles.ASSL.XOffset = get(hObject, 'Max');
end

axes(handles.ASSLSpecAxis);
handles.ASSL.SpecAxisLimits = [handles.ASSL.XOffset (handles.ASSL.XOffset + (handles.ASSL.SpecAxisLimits(2) - handles.ASSL.SpecAxisLimits(1))) handles.ASSL.SpecAxisLimits(3:4)];
axis(handles.ASSL.SpecAxisLimits);

axes(handles.ASSLAmpAxis);
handles.ASSL.AmpAxisLimits = [handles.ASSL.XOffset (handles.ASSL.XOffset + (handles.ASSL.AmpAxisLimits(2) - handles.ASSL.AmpAxisLimits(1))) handles.ASSL.AmpAxisLimits(3:4)];
axis(handles.ASSL.AmpAxisLimits);

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function XOffsetSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to XOffsetSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in SingleMultipleThresholdToggle.
function SingleMultipleThresholdToggle_Callback(hObject, eventdata, handles)
% hObject    handle to SingleMultipleThresholdToggle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of SingleMultipleThresholdToggle


% --- Executes on button press in ManualThresholdButton.
function ManualThresholdButton_Callback(hObject, eventdata, handles)
% hObject    handle to ManualThresholdButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ManualThresholdButton
if (get(hObject, 'Value') == 1)
    set(handles.ThresholdEdit, 'Enable', 'on');
    handles.ASSL.AutoThreshold = 0;
    set(handles.AutoThresholdButton, 'Value', 0);
else
    set(handles.ThresholdEdit, 'Enable', 'off');
    handles.ASSL.AutoThreshold = 1;
    set(handles.AutoThresholdButton, 'Value', 1);    
end
guidata(hObject, handles);

% --- Executes on button press in AutoThresholdButton.
function AutoThresholdButton_Callback(hObject, eventdata, handles)
% hObject    handle to AutoThresholdButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of AutoThresholdButton
if (get(hObject, 'Value') == 1)
    set(handles.ThresholdEdit, 'Enable', 'off');
    handles.ASSL.AutoThreshold = 0;
    set(handles.ManualThresholdButton, 'Value', 0);
else
    set(handles.ThresholdEdit, 'Enable', 'on');
    handles.ASSL.AutoThreshold = 0;
    set(handles.ManualThresholdButton, 'Value', 1);    
end
guidata(hObject, handles);


function ThresholdEdit_Callback(hObject, eventdata, handles)
% hObject    handle to ThresholdEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ThresholdEdit as text
%        str2double(get(hObject,'String')) returns contents of ThresholdEdit as a double
handles.ASSL.FixedThreshold = str2double(get(hObject, 'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function ThresholdEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ThresholdEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in ToBeUsedFeaturesListBox.
function ToBeUsedFeaturesListBox_Callback(hObject, eventdata, handles)
% hObject    handle to ToBeUsedFeaturesListBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ToBeUsedFeaturesListBox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ToBeUsedFeaturesListBox


% --- Executes during object creation, after setting all properties.
function ToBeUsedFeaturesListBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ToBeUsedFeaturesListBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function FFTWinSizeTempMatchEdit_Callback(hObject, eventdata, handles)
% hObject    handle to FFTWinSizeTempMatchEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FFTWinSizeTempMatchEdit as text
%        str2double(get(hObject,'String')) returns contents of FFTWinSizeTempMatchEdit as a double
handles.ASSL.FFTWinSizeTempMatch = str2double(get(hObject, 'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function FFTWinSizeTempMatchEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FFTWinSizeTempMatchEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in SetTemplateDirButton.
function SetTemplateDirButton_Callback(hObject, eventdata, handles)
% hObject    handle to SetTemplateDirButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.ASSL.TemplateDir = uigetdir(pwd, 'Choose the directory with the template files');

guidata(hObject, handles);


% --- Executes on button press in ReviewLabelsPushButton.
function ReviewLabelsPushButton_Callback(hObject, eventdata, handles)
% hObject    handle to ReviewLabelsPushButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ASSLReviewTemplateMatching(handles.ASSL);
guidata(hObject, handles);


% --- Executes on button press in FixSyllBoundariesButton.
function FixSyllBoundariesButton_Callback(hObject, eventdata, handles)
% hObject    handle to FixSyllBoundariesButton (see GCBO)
% eventdata  reserved - to be defined in a future version of M
% handles    empty - handles not created until after all CreateFcns called
% handles    empty - handles not created until after all CreateFcns calledATLAB
% handles    structure with handles and user data (see GUIDATA)

ASSLAdjustSyllableBoundaries(handles);

guidata(hObject, handles);


% --- Executes on button press in SegementSongsAronovFeeButton.
function SegementSongsAronovFeeButton_Callback(hObject, eventdata, handles)
% hObject    handle to SegementSongsAronovFeeButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

for i = 1:length(handles.ASSL.FileName),
    if (isfield(handles.ASSL, 'SyllLabels'))
        if (length(handles.ASSL.SyllLabels) >= i)
            if (~isempty(handles.ASSL.SyllLabels{i}))
                handles.ASSL.SyllLabels(i) = [];
            end
        end
    end
    
    if (exist([handles.ASSL.NoteFileDirName, handles.ASSL.FileName{i}, '.not.mat'], 'file'))
        Temp = load([handles.ASSL.NoteFileDirName, handles.ASSL.FileName{i}, '.not.mat']);
        handles.ASSL.SyllOnsets{i} = Temp.onsets;
        handles.ASSL.SyllOffsets{i} = Temp.offsets;
        handles.ASSL.SyllLabels{i} = Temp.labels;
        if (isfield(Temp, 'threshold'))
            handles.ASSL.Threshold{i} = Temp.threshold;
        else
            handles.ASSL.Threshold{i} = -30;
        end
        if (~isfield(Temp, 'FileLength'))
            [RawData, Fs] = ASSLGetRawData(handles.ASSL.DirName, handles.ASSL.FileName{i}, handles.ASSL.FileType, handles.ASSL.SongChanNo);
            handles.ASSL.FileDur{i} = length(RawData)/Fs;
        end
        disp(['Loaded syllable data from ', handles.ASSL.FileName{i}, '.not.mat']);
        continue;
    end
    [RawData, Fs] = ASSLGetRawData(handles.ASSL.DirName, handles.ASSL.FileName{i}, handles.ASSL.FileType, handles.ASSL.SongChanNo);
    handles.ASSL.FileDur{i} = length(RawData)/Fs;
    
    Time = (1:1:length(RawData))/Fs;
    
    [LogAmplitude] = ASSLCalculateLogAmplitudeAronovFee(RawData, Fs, Time, handles.ASSL.FFTWinSizeSegmenting, handles.ASSL.FFTWinOverlapSegmenting);
    
    if (handles.ASSL.AutoThreshold == 0)
        handles.ASSL.Threshold{i} = handles.ASSL.FixedThreshold;
    else
        %Obj = gmdistribution.fit(LogAmplitude', 2);
        handles.ASSL.Threshold{i} = ASSLCalculateFisherThreshold(LogAmplitude);
        %handles.ASSL.Threshold{i} = mean(Obj.mu);
    end
    
    [handles.ASSL.SyllOnsets{i}, handles.ASSL.SyllOffsets{i}] = ASSLSegmentDataAronovFee(LogAmplitude, Fs, handles.ASSL.MinInt, handles.ASSL.MinDur, handles.ASSL.Threshold{i});
    for j = 1:length(handles.ASSL.SyllOnsets{i}),
        handles.ASSL.SyllLabels{i}(j) = '0';
    end
    disp([' Detected ', num2str(length(handles.ASSL.SyllOnsets{i})), ' syllables for ', handles.ASSL.FileName{i}]);
end

set(handles.SongFileNameText, 'String', ['Song File Name : ', handles.ASSL.FileName{handles.ASSL.FileIndex}]);
[RawData, Fs] = ASSLGetRawData(handles.ASSL.DirName, handles.ASSL.FileName{handles.ASSL.FileIndex}, handles.ASSL.FileType, handles.ASSL.SongChanNo);
    
Time = (1:1:length(RawData))/Fs;
[LogAmplitude] = ASSLCalculateLogAmplitudeAronovFee(RawData, Fs, Time, handles.ASSL.FFTWinSizeSegmenting, handles.ASSL.FFTWinOverlapSegmenting);

[handles.ASSL.SpecAxisLimits, handles.ASSL.AmpAxisLimits] = ASSLPlotData(handles.ASSL.DirName, handles.ASSL.FileName{handles.ASSL.FileIndex}, handles.ASSL.FileType, Time, LogAmplitude, handles.ASSLSpecAxis, handles.ASSLAmpAxis, handles.ASSL.Threshold{handles.ASSL.FileIndex}, handles.ASSL.SyllOnsets{handles.ASSL.FileIndex}, handles.ASSL.SyllOffsets{handles.ASSL.FileIndex});

set(handles.XOffsetSlider, 'Max', Time(end));
set(handles.XZoomSlider, 'Max', Time(end));
set(handles.XZoomSlider, 'Min', 0.01);
set(handles.XZoomSlider, 'Value', get(handles.XZoomSlider, 'Max'));
    
guidata(hObject, handles);


% --- Executes on button press in LabelSongsButton.
function LabelSongsButton_Callback(hObject, eventdata, handles)
% hObject    handle to LabelSongsButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ASSLLabelSongs(handles.ASSL);
guidata(hObject, handles);

% --- Executes on button press in AnalyzeSyllablesButton.
function AnalyzeSyllablesButton_Callback(hObject, eventdata, handles)
% hObject    handle to AnalyzeSyllablesButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ASSLAnalyzeSyllables(handles.ASSL);
guidata(hObject, handles);


% --- Executes on button press in SaveDataButton.
function SaveDataButton_Callback(hObject, eventdata, handles)
% hObject    handle to SaveDataButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Filesep = filesep;

for i = 1:length(handles.ASSL.ToBeUsedFeatures),
    if (handles.ASSL.NoteFileDirName(end) == Filesep)
        if ((isfield(handles.ASSL, 'FileListName')) && (~isempty(handles.ASSL.FileListName)))
            RootOutputFileName = [handles.ASSL.NoteFileDirName, handles.ASSL.FileListName];
        else
            RootOutputFileName = [handles.ASSL.NoteFileDirName, handles.ASSL.FileName{1}];
        end
    else
        if ((isfield(handles.ASSL, 'FileListName')) && (~isempty(handles.ASSL.FileListName)))
            RootOutputFileName = [handles.ASSL.NoteFileDirName, Filesep, handles.ASSL.FileListName];
        else
            RootOutputFileName = [handles.ASSL.NoteFileDirName, Filesep, handles.ASSL.FileName{1}];
        end
    end
end

OutputFileName = [RootOutputFileName, '.ASSLData.mat'];
save(OutputFileName, 'handles');
disp('Finished saving data');
guidata(hObject, handles);


% --- Executes on button press in SaveExcelButton.
function SaveExcelButton_Callback(hObject, eventdata, handles)
% hObject    handle to SaveExcelButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Filesep = filesep;

if (handles.ASSL.NoteFileDirName(end) == Filesep)
    if ((isfield(handles.ASSL, 'FileListName')) && (~isempty(handles.ASSL.FileListName)))
        RootOutputFileName = [handles.ASSL.NoteFileDirName, handles.ASSL.FileListName];
    else
        RootOutputFileName = [handles.ASSL.NoteFileDirName, handles.ASSL.FileName{1}];
    end
else
    if ((isfield(handles.ASSL, 'FileListName')) && (~isempty(handles.ASSL.FileListName)))
        RootOutputFileName = [handles.ASSL.NoteFileDirName, Filesep, handles.ASSL.FileListName];
    else
        RootOutputFileName = [handles.ASSL.NoteFileDirName, Filesep, handles.ASSL.FileName{1}];
    end
end

TextOutputFileName = [RootOutputFileName, '.ASSLOnsetOffsetData.txt'];
OutputFileName = [RootOutputFileName, '.ASSLOnsetOffsetData.xls'];

Fid = fopen(TextOutputFileName, 'w');
RowIndex = 1;
for i = 1:length(handles.ASSL.SyllOnsets),
    for j = 1:length(handles.ASSL.SyllOnsets{i}),
        fprintf(Fid, '%i\t%c\t%g\t%g\n', i, handles.ASSL.SyllLabels{i}(j), handles.ASSL.SyllOnsets{i}(j), handles.ASSL.SyllOffsets{i}(j));
        Temp(RowIndex,:) = {i, handles.ASSL.SyllLabels{i}(j), handles.ASSL.SyllOnsets{i}(j), handles.ASSL.SyllOffsets{i}(j)};
        RowIndex = RowIndex + 1;
    end
end
fclose(Fid);
disp(['Wrote data about labels, onsets and offsets to ', TextOutputFileName]); 
xlswrite(OutputFileName, Temp, 1, 'A1');
disp(['Wrote data about labels, onsets and offsets to ', OutputFileName]); 


% --- Executes on button press in DelFeatValueFilesButton.
function DelFeatValueFilesButton_Callback(hObject, eventdata, handles)
% hObject    handle to DelFeatValueFilesButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Filesep = filesep;

for i = 1:length(handles.ASSL.ToBeUsedFeatures),
    if (handles.ASSL.NoteFileDirName(end) == Filesep)
        if ((isfield(handles.ASSL, 'FileListName')) && (~isempty(handles.ASSL.FileListName)))
            RootOutputFileName = [handles.ASSL.NoteFileDirName, handles.ASSL.FileListName];
        else
            RootOutputFileName = [handles.ASSL.NoteFileDirName, handles.ASSL.FileName{1}];
        end
    else
        if ((isfield(handles.ASSL, 'FileListName')) && (~isempty(handles.ASSL.FileListName)))
            RootOutputFileName = [handles.ASSL.NoteFileDirName, Filesep, handles.ASSL.FileListName];
        else
            RootOutputFileName = [handles.ASSL.NoteFileDirName, Filesep, handles.ASSL.FileName{1}];
        end
    end
    
    FeatValuesFileName = [RootOutputFileName, '.FeatValues.mat'];
    RawFeatValuesFileName = [RootOutputFileName, '.RawFeatValues.mat'];
    SyllIndicesFileName = [RootOutputFileName, '.SyllIndices.mat'];
    SyllIndexLabelsFileName = [RootOutputFileName, '.SyllIndexLabels.mat'];
    RawOutputFileName = [RootOutputFileName, '.', handles.ASSL.ToBeUsedFeatures{i}, '.raw.mat'];
    OutputFileName = [RootOutputFileName, '.', handles.ASSL.ToBeUsedFeatures{i}, '.mat'];
    
    delete(OutputFileName, FeatValuesFileName, RawFeatValuesFileName, SyllIndicesFileName, SyllIndexLabelsFileName, RawOutputFileName);
end

disp('Finished deleting feature value files');


% --- Executes on button press in ChooseSyllFFBoundariesButton.
function ChooseSyllFFBoundariesButton_Callback(hObject, eventdata, handles)
% hObject    handle to ChooseSyllFFBoundariesButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ASSLChooseSyllFFBoundaries(handles.ASSL);
guidata(hObject, handles);
