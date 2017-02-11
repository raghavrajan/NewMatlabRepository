function varargout = ASSLPickBouts(varargin)
%ASSLPICKBOUTSMAINFIG M-file for ASSLPickBoutsMainFig.fig
%      ASSLPICKBOUTSMAINFIG, by itself, creates a new ASSLPICKBOUTSMAINFIG or raises the existing
%      singleton*.
%
%      H = ASSLPICKBOUTSMAINFIG returns the handle to a new ASSLPICKBOUTSMAINFIG or the handle to
%      the existing singleton*.
%
%      ASSLPICKBOUTSMAINFIG('Property','Value',...) creates a new ASSLPICKBOUTSMAINFIG using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to ASSLPickBouts_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      ASSLPICKBOUTSMAINFIG('CALLBACK') and ASSLPICKBOUTSMAINFIG('CALLBACK',hObject,...) call the
%      local function named CALLBACK in ASSLPICKBOUTSMAINFIG.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ASSLPickBoutsMainFig

% Last Modified by GUIDE v2.5 23-Apr-2016 18:02:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ASSLPickBouts_OpeningFcn, ...
                   'gui_OutputFcn',  @ASSLPickBouts_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
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


% --- Executes just before ASSLPickBoutsMainFig is made visible.
function ASSLPickBouts_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

handles.ASSLPickBoutsMainFig.TimeStep = str2double(get(handles.TimeStepEdit, 'String'));
handles.ASSLPickBoutsMainFig.InterBoutInterval = str2double(get(handles.InterBoutIntervalEdit, 'String'))*1000;
handles.ASSLPickBoutsMainFig.BoutPaddingTime = str2double(get(handles.BoutPaddingTimeEdit, 'String'))*1000;

handles.ASSLPickBoutsMainFig.LoResHiRes = 1;
set(handles.LoResHiResToggle, 'String', 'High Res Spectrogram');
set(handles.LoResHiResToggle, 'Value', handles.ASSLPickBoutsMainFig.LoResHiRes);

if (nargin >= 1)
    handles.ASSLPickBoutsMainFig = varargin{1};
    handles.ASSLPickBoutsMainFig.LoResHiRes = 1;
    handles.ASSLPickBoutsMainFig.TimeStep = str2double(get(handles.TimeStepEdit, 'String'));
    handles.ASSLPickBoutsMainFig.InterBoutInterval = str2double(get(handles.InterBoutIntervalEdit, 'String'))*1000;
    handles.ASSLPickBoutsMainFig.BoutPaddingTime = str2double(get(handles.BoutPaddingTimeEdit, 'String'))*1000;
    
    for i = 1:length(handles.ASSLPickBoutsMainFig.FileName),
        [handles.BoutOnsetsOffsets{i}] = ASSLGetBoutEdges(handles.ASSLPickBoutsMainFig.SyllOnsets{i}, handles.ASSLPickBoutsMainFig.SyllOffsets{i}, handles.ASSLPickBoutsMainFig.InterBoutInterval, handles.ASSLPickBoutsMainFig.BoutPaddingTime, handles.ASSLPickBoutsMainFig.FileDur{i});
    end
    
    handles.ASSLPickBoutsMainFig.FileIndex = 1;
    set(handles.SongFileNameTextLabel, 'String', ['Song File Name : ', handles.ASSLPickBoutsMainFig.FileName{handles.ASSLPickBoutsMainFig.FileIndex}, ' : #', num2str(handles.ASSLPickBoutsMainFig.FileIndex), ' of ', num2str(length(handles.ASSLPickBoutsMainFig.FileName)), ' files']);
    
    [RawData, Fs] = ASSLGetRawData(handles.ASSLPickBoutsMainFig.DirName, handles.ASSLPickBoutsMainFig.FileName{handles.ASSLPickBoutsMainFig.FileIndex}, handles.ASSLPickBoutsMainFig.FileType, handles.ASSLPickBoutsMainFig.SongChanNo);

    Time = (1:1:length(RawData))/Fs;

    [LogAmplitude] = ASSLCalculateLogAmplitude(RawData, Fs, Time, handles.ASSLPickBoutsMainFig.FFTWinSizeSegmenting, handles.ASSLPickBoutsMainFig.FFTWinOverlapSegmenting);

%    set(handles.SongFileNameText, 'String', ['Song File Name : ', handles.ASSLPickBoutsMainFig.FileName{handles.ASSLPickBoutsMainFig.FileIndex}]);

    [handles.ASSLPickBoutsMainFig.SpecAxisLimits, handles.ASSLPickBoutsMainFig.LabelAxisLimits, handles.ASSLPickBoutsMainFig.AmpAxisLimits] = ASSLReviewTMPlotData(handles, Time, LogAmplitude);
    
    if (length(handles.ASSLPickBoutsMainFig.SyllOnsets{handles.ASSLPickBoutsMainFig.FileIndex}) > 1)
        axes(handles.ReviewSpecAxis);
        handles.ASSLPickBoutsMainFig.BoutPlotLine = plot([handles.BoutOnsetsOffsets{handles.ASSLPickBoutsMainFig.FileIndex}]'/1000, ones(size([handles.BoutOnsetsOffsets{handles.ASSLPickBoutsMainFig.FileIndex}]))'*7500, 'b', 'LineWidth', 3);
    end
end

handles.ASSLPickBoutsMainFig.ZoomSpecAxisLimits = handles.ASSLPickBoutsMainFig.SpecAxisLimits;
handles.ASSLPickBoutsMainFig.ZoomAmpAxisLimits = handles.ASSLPickBoutsMainFig.AmpAxisLimits;
handles.ASSLPickBoutsMainFig.ZoomLabelAxisLimits = handles.ASSLPickBoutsMainFig.LabelAxisLimits;

% Update handles structure
guidata(hObject, handles);

UpdateBoutNumber(handles);

% UIWAIT makes ASSLPickBoutsMainFig wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ASSLPickBouts_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in SaveBoutNoteFilesButton.
function SaveBoutNoteFilesButton_Callback(hObject, eventdata, handles)
% hObject    handle to SaveBoutNoteFilesButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

FileSep = filesep;
DirName = uigetdir(pwd, 'Pick the directory where bout files should be saved as .not.mat files in BoutNoteFiles dir');
PresentDir = pwd;
cd(DirName);
if (~exist('BoutNoteFiles', 'dir'))
    mkdir('BoutNoteFiles');
end
cd('BoutNoteFiles');
for i = 1:length(handles.ASSLPickBoutsMainFig.FileName),
    if (~isempty(handles.BoutOnsetsOffsets{i}))
        onsets = handles.BoutOnsetsOffsets{i}(:,1);
        offsets = handles.BoutOnsetsOffsets{i}(:,2);
        labels = repmat('B', 1, length(onsets));
    else
        onsets = [];
        offsets = [];
        labels = [];
    end
    min_int = handles.ASSLPickBoutsMainFig.InterBoutInterval;
    save([handles.ASSLPickBoutsMainFig.FileName{i}, '.not.mat'], 'onsets', 'offsets', 'labels', 'min_int');
end
cd(PresentDir);
msgbox(['Finished writing bout .not.mat files to ', DirName, FileSep, 'BoutNoteFiles']);
guidata(hObject, handles);

function TimeStepEdit_Callback(hObject, eventdata, handles)
% hObject    handle to TimeStepEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TimeStepEdit as text
%        str2double(get(hObject,'String')) returns contents of TimeStepEdit as a double
handles.ASSLPickBoutsMainFig.TimeStep = str2double(get(hObject, 'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function TimeStepEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TimeStepEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in PickBoutEdgesButton.
function PickBoutEdgesButton_Callback(hObject, eventdata, handles)
% hObject    handle to PickBoutEdgesButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Flag = 1;

set(handles.InstructionsTextLabel, 'String', 'Instructions: Bouts are marked with blue lines on top. Left click to delete a bout. Be sure to click within the syllable. Type N for next time step, P for previous time step, n for next file and p for previous file. Type q to quit');
axes(handles.ReviewSpecAxis);

while (Flag)
    [x, y, button] = ginput(1);
    switch button
        case 113 % typing q
            break;
            
        case 1 % left click
            x(1) = x(1) * 1000;
            TempBouts = handles.BoutOnsetsOffsets{handles.ASSLPickBoutsMainFig.FileIndex};
            if (~isempty(TempBouts)),
                for i = 1:size(TempBouts,1),
                    if ((x(1) >= TempBouts(i,1)) && (x(1) <= TempBouts(i,2)))
                        TempBouts(i,:) = [];
                        break;
                    end
                end
                handles.BoutOnsetsOffsets{handles.ASSLPickBoutsMainFig.FileIndex} = TempBouts;
                axes(handles.ReviewSpecAxis);
                delete(handles.ASSLPickBoutsMainFig.BoutPlotLine);
                handles.ASSLPickBoutsMainFig.BoutPlotLine = plot([handles.BoutOnsetsOffsets{handles.ASSLPickBoutsMainFig.FileIndex}]'/1000, ones(size([handles.BoutOnsetsOffsets{handles.ASSLPickBoutsMainFig.FileIndex}]))'*7500, 'b', 'LineWidth', 3);
                UpdateBoutNumber(handles);
            end
        
        case 78 % typing N
            NextTimeButton_Callback(hObject, eventdata, handles);
            handles.ASSLPickBoutsMainFig = getappdata(findobj('Tag', 'ASSLPickBouts'), 'Data');

        case 80 % typing P
            PrevTimeButton_Callback(hObject, eventdata, handles);
            handles.ASSLPickBoutsMainFig = getappdata(findobj('Tag', 'ASSLPickBouts'), 'Data');
            
        case 110 % typing n
            NextFileButton_Callback(hObject, eventdata, handles);
            handles.ASSLPickBoutsMainFig = getappdata(findobj('Tag', 'ASSLPickBouts'), 'Data');
            
        case 112 % typing p
            PrevFileButton_Callback(hObject, eventdata, handles);
            handles.ASSLPickBoutsMainFig = getappdata(findobj('Tag', 'ASSLPickBouts'), 'Data');
            
    end
end
guidata(hObject, handles);


% --- Executes on button press in NextFileButton.
function NextFileButton_Callback(hObject, eventdata, handles)
% hObject    handle to NextFileButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.ASSLPickBoutsMainFig.FileIndex = handles.ASSLPickBoutsMainFig.FileIndex + 1;
if (handles.ASSLPickBoutsMainFig.FileIndex < 1)
    handles.ASSLPickBoutsMainFig.FileIndex = 1;
end

if (handles.ASSLPickBoutsMainFig.FileIndex > length(handles.ASSLPickBoutsMainFig.FileName))
    handles.ASSLPickBoutsMainFig.FileIndex = length(handles.ASSLPickBoutsMainFig.FileName);
end

set(handles.SongFileNameTextLabel, 'String', ['Song File Name : ', handles.ASSLPickBoutsMainFig.FileName{handles.ASSLPickBoutsMainFig.FileIndex}, ' : #', num2str(handles.ASSLPickBoutsMainFig.FileIndex), ' of ', num2str(length(handles.ASSLPickBoutsMainFig.FileName)), ' files']);
[RawData, Fs] = ASSLGetRawData(handles.ASSLPickBoutsMainFig.DirName, handles.ASSLPickBoutsMainFig.FileName{handles.ASSLPickBoutsMainFig.FileIndex}, handles.ASSLPickBoutsMainFig.FileType, handles.ASSLPickBoutsMainFig.SongChanNo);
Time = (1:1:length(RawData))/Fs;
[LogAmplitude] = ASSLCalculateLogAmplitude(RawData, Fs, Time, handles.ASSLPickBoutsMainFig.FFTWinSizeSegmenting, handles.ASSLPickBoutsMainFig.FFTWinOverlapSegmenting);

[handles.ASSLPickBoutsMainFig.SpecAxisLimits, handles.ASSLPickBoutsMainFig.LabelAxisLimits, handles.ASSLPickBoutsMainFig.AmpAxisLimits] = ASSLReviewTMPlotData(handles, Time, LogAmplitude);

if (length(handles.ASSLPickBoutsMainFig.SyllOnsets{handles.ASSLPickBoutsMainFig.FileIndex}) > 1)
    axes(handles.ReviewSpecAxis);
    handles.ASSLPickBoutsMainFig.BoutPlotLine = plot([handles.BoutOnsetsOffsets{handles.ASSLPickBoutsMainFig.FileIndex}]'/1000, ones(size([handles.BoutOnsetsOffsets{handles.ASSLPickBoutsMainFig.FileIndex}]))'*7500, 'b', 'LineWidth', 3);
end

handles.ASSLPickBoutsMainFig.ZoomSpecAxisLimits = handles.ASSLPickBoutsMainFig.SpecAxisLimits;
handles.ASSLPickBoutsMainFig.ZoomAmpAxisLimits = handles.ASSLPickBoutsMainFig.AmpAxisLimits;
handles.ASSLPickBoutsMainFig.ZoomLabelAxisLimits = handles.ASSLPickBoutsMainFig.LabelAxisLimits;

setappdata(findobj('Tag', 'ASSLPickBouts'), 'Data', handles.ASSLPickBoutsMainFig);
guidata(hObject, handles);

UpdateBoutNumber(handles);

% --- Executes on button press in PrevFileButton.
function PrevFileButton_Callback(hObject, eventdata, handles)
% hObject    handle to PrevFileButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.ASSLPickBoutsMainFig.FileIndex = handles.ASSLPickBoutsMainFig.FileIndex - 1;
if (handles.ASSLPickBoutsMainFig.FileIndex < 1)
    handles.ASSLPickBoutsMainFig.FileIndex = 1;
end

if (handles.ASSLPickBoutsMainFig.FileIndex > length(handles.ASSLPickBoutsMainFig.FileName))
    handles.ASSLPickBoutsMainFig.FileIndex = length(handles.ASSLPickBoutsMainFig.FileName);
end

set(handles.SongFileNameTextLabel, 'String', ['Song File Name : ', handles.ASSLPickBoutsMainFig.FileName{handles.ASSLPickBoutsMainFig.FileIndex}, ' : #', num2str(handles.ASSLPickBoutsMainFig.FileIndex), ' of ', num2str(length(handles.ASSLPickBoutsMainFig.FileName)), ' files']);
[RawData, Fs] = ASSLGetRawData(handles.ASSLPickBoutsMainFig.DirName, handles.ASSLPickBoutsMainFig.FileName{handles.ASSLPickBoutsMainFig.FileIndex}, handles.ASSLPickBoutsMainFig.FileType, handles.ASSLPickBoutsMainFig.SongChanNo);
Time = (1:1:length(RawData))/Fs;
[LogAmplitude] = ASSLCalculateLogAmplitude(RawData, Fs, Time, handles.ASSLPickBoutsMainFig.FFTWinSizeSegmenting, handles.ASSLPickBoutsMainFig.FFTWinOverlapSegmenting);

[handles.ASSLPickBoutsMainFig.SpecAxisLimits, handles.ASSLPickBoutsMainFig.LabelAxisLimits, handles.ASSLPickBoutsMainFig.AmpAxisLimits] = ASSLReviewTMPlotData(handles, Time, LogAmplitude);

if (length(handles.ASSLPickBoutsMainFig.SyllOnsets{handles.ASSLPickBoutsMainFig.FileIndex}) > 1)
    axes(handles.ReviewSpecAxis);
    handles.ASSLPickBoutsMainFig.BoutPlotLine = plot([handles.BoutOnsetsOffsets{handles.ASSLPickBoutsMainFig.FileIndex}]'/1000, ones(size([handles.BoutOnsetsOffsets{handles.ASSLPickBoutsMainFig.FileIndex}]))'*7500, 'b', 'LineWidth', 3);
end

handles.ASSLPickBoutsMainFig.ZoomSpecAxisLimits = handles.ASSLPickBoutsMainFig.SpecAxisLimits;
handles.ASSLPickBoutsMainFig.ZoomAmpAxisLimits = handles.ASSLPickBoutsMainFig.AmpAxisLimits;
handles.ASSLPickBoutsMainFig.ZoomLabelAxisLimits = handles.ASSLPickBoutsMainFig.LabelAxisLimits;

setappdata(findobj('Tag', 'ASSLPickBouts'), 'Data', handles.ASSLPickBoutsMainFig);
guidata(hObject, handles);

UpdateBoutNumber(handles);

% --- Executes on button press in NextTimeButton.
function NextTimeButton_Callback(hObject, eventdata, handles)
% hObject    handle to NextTimeButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ((handles.ASSLPickBoutsMainFig.ZoomSpecAxisLimits(2) - handles.ASSLPickBoutsMainFig.ZoomSpecAxisLimits(1)) > 1.1*handles.ASSLPickBoutsMainFig.TimeStep)
    handles.ASSLPickBoutsMainFig.ZoomSpecAxisLimits(2) = handles.ASSLPickBoutsMainFig.SpecAxisLimits(1) + handles.ASSLPickBoutsMainFig.TimeStep;
else
    handles.ASSLPickBoutsMainFig.ZoomSpecAxisLimits(1) = handles.ASSLPickBoutsMainFig.ZoomSpecAxisLimits(1) + 0.9*handles.ASSLPickBoutsMainFig.TimeStep;
    handles.ASSLPickBoutsMainFig.ZoomSpecAxisLimits(2) = handles.ASSLPickBoutsMainFig.ZoomSpecAxisLimits(1) + handles.ASSLPickBoutsMainFig.TimeStep;
end

if (handles.ASSLPickBoutsMainFig.ZoomSpecAxisLimits(1) < 0)
    handles.ASSLPickBoutsMainFig.ZoomSpecAxisLimits(1) = 0;
    handles.ASSLPickBoutsMainFig.ZoomSpecAxisLimits(2) = handles.ASSLPickBoutsMainFig.ZoomSpecAxisLimits(1) + handles.ASSLPickBoutsMainFig.TimeStep;
end

if (handles.ASSLPickBoutsMainFig.ZoomSpecAxisLimits(2) > handles.ASSLPickBoutsMainFig.SpecAxisLimits(2))
    handles.ASSLPickBoutsMainFig.ZoomSpecAxisLimits(2) = handles.ASSLPickBoutsMainFig.SpecAxisLimits(2) + handles.ASSLPickBoutsMainFig.TimeStep*0.1;
    handles.ASSLPickBoutsMainFig.ZoomSpecAxisLimits(1) = handles.ASSLPickBoutsMainFig.ZoomSpecAxisLimits(2) - handles.ASSLPickBoutsMainFig.TimeStep;
end

handles.ASSLPickBoutsMainFig.ZoomAmpAxisLimits(1:2) = handles.ASSLPickBoutsMainFig.ZoomSpecAxisLimits(1:2);
handles.ASSLPickBoutsMainFig.ZoomLabelAxisLimits(1:2) = handles.ASSLPickBoutsMainFig.ZoomSpecAxisLimits(1:2);

axes(handles.ReviewSpecAxis);
axis(handles.ASSLPickBoutsMainFig.ZoomSpecAxisLimits);

axes(handles.ReviewLabelAxis);
axis(handles.ASSLPickBoutsMainFig.ZoomLabelAxisLimits);

axes(handles.ReviewAmplitudeAxis);
axis(handles.ASSLPickBoutsMainFig.ZoomAmpAxisLimits);

setappdata(findobj('Tag', 'ASSLPickBouts'), 'Data', handles.ASSLPickBoutsMainFig);
guidata(hObject, handles);

% --- Executes on button press in PrevTimeButton.
function PrevTimeButton_Callback(hObject, eventdata, handles)
% hObject    handle to PrevTimeButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


if ((handles.ASSLPickBoutsMainFig.ZoomSpecAxisLimits(2) - handles.ASSLPickBoutsMainFig.ZoomSpecAxisLimits(1)) > 1.1*handles.ASSLPickBoutsMainFig.TimeStep)
    handles.ASSLPickBoutsMainFig.ZoomSpecAxisLimits(2) = handles.ASSLPickBoutsMainFig.SpecAxisLimits(1) + handles.ASSLPickBoutsMainFig.TimeStep;
else
    handles.ASSLPickBoutsMainFig.ZoomSpecAxisLimits(1) = handles.ASSLPickBoutsMainFig.ZoomSpecAxisLimits(1) - 0.9*handles.ASSLPickBoutsMainFig.TimeStep;
    handles.ASSLPickBoutsMainFig.ZoomSpecAxisLimits(2) = handles.ASSLPickBoutsMainFig.ZoomSpecAxisLimits(1) + handles.ASSLPickBoutsMainFig.TimeStep;
end

if (handles.ASSLPickBoutsMainFig.ZoomSpecAxisLimits(1) < 0)
    handles.ASSLPickBoutsMainFig.ZoomSpecAxisLimits(1) = 0;
    handles.ASSLPickBoutsMainFig.ZoomSpecAxisLimits(2) = handles.ASSLPickBoutsMainFig.ZoomSpecAxisLimits(1) + handles.ASSLPickBoutsMainFig.TimeStep;
end

if (handles.ASSLPickBoutsMainFig.ZoomSpecAxisLimits(2) > handles.ASSLPickBoutsMainFig.SpecAxisLimits(2))
    handles.ASSLPickBoutsMainFig.ZoomSpecAxisLimits(2) = handles.ASSLPickBoutsMainFig.SpecAxisLimits(2) + handles.ASSLPickBoutsMainFig.TimeStep*0.1;
    handles.ASSLPickBoutsMainFig.ZoomSpecAxisLimits(1) = handles.ASSLPickBoutsMainFig.ZoomSpecAxisLimits(2) - handles.ASSLPickBoutsMainFig.TimeStep;
end


handles.ASSLPickBoutsMainFig.ZoomAmpAxisLimits(1:2) = handles.ASSLPickBoutsMainFig.ZoomSpecAxisLimits(1:2);
handles.ASSLPickBoutsMainFig.ZoomLabelAxisLimits(1:2) = handles.ASSLPickBoutsMainFig.ZoomSpecAxisLimits(1:2);

axes(handles.ReviewSpecAxis);
axis(handles.ASSLPickBoutsMainFig.ZoomSpecAxisLimits);

axes(handles.ReviewLabelAxis);
axis(handles.ASSLPickBoutsMainFig.ZoomLabelAxisLimits);

axes(handles.ReviewAmplitudeAxis);
axis(handles.ASSLPickBoutsMainFig.ZoomAmpAxisLimits);

setappdata(findobj('Tag', 'ASSLPickBouts'), 'Data', handles.ASSLPickBoutsMainFig);
guidata(hObject, handles);


% --- Executes on button press in JumpFileButton.
function JumpFileButton_Callback(hObject, eventdata, handles)
% hObject    handle to JumpFileButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

NewFileIndex = inputdlg('Enter the file # to jump to', 'File #');
NewFileIndex = str2double(NewFileIndex{1});

handles.ASSLPickBoutsMainFig.FileIndex = NewFileIndex;
if (handles.ASSLPickBoutsMainFig.FileIndex < 1)
    handles.ASSLPickBoutsMainFig.FileIndex = 1;
end

if (handles.ASSLPickBoutsMainFig.FileIndex > length(handles.ASSLPickBoutsMainFig.FileName))
    handles.ASSLPickBoutsMainFig.FileIndex = length(handles.ASSLPickBoutsMainFig.FileName);
end

set(handles.SongFileNameTextLabel, 'String', ['Song File Name : ', handles.ASSLPickBoutsMainFig.FileName{handles.ASSLPickBoutsMainFig.FileIndex}, ' : #', num2str(handles.ASSLPickBoutsMainFig.FileIndex), ' of ', num2str(length(handles.ASSLPickBoutsMainFig.FileName)), ' files']);
[RawData, Fs] = ASSLGetRawData(handles.ASSLPickBoutsMainFig.DirName, handles.ASSLPickBoutsMainFig.FileName{handles.ASSLPickBoutsMainFig.FileIndex}, handles.ASSLPickBoutsMainFig.FileType, handles.ASSLPickBoutsMainFig.SongChanNo);
Time = (1:1:length(RawData))/Fs;
[LogAmplitude] = ASSLCalculateLogAmplitude(RawData, Fs, Time, handles.ASSLPickBoutsMainFig.FFTWinSizeSegmenting, handles.ASSLPickBoutsMainFig.FFTWinOverlapSegmenting);

[handles.ASSLPickBoutsMainFig.SpecAxisLimits, handles.ASSLPickBoutsMainFig.LabelAxisLimits, handles.ASSLPickBoutsMainFig.AmpAxisLimits] = ASSLReviewTMPlotData(handles, Time, LogAmplitude);

if (length(handles.ASSLPickBoutsMainFig.SyllOnsets{handles.ASSLPickBoutsMainFig.FileIndex}) > 1)
    axes(handles.ReviewSpecAxis);
    handles.ASSLPickBoutsMainFig.BoutPlotLine = plot([handles.BoutOnsetsOffsets{handles.ASSLPickBoutsMainFig.FileIndex}]'/1000, ones(size([handles.BoutOnsetsOffsets{handles.ASSLPickBoutsMainFig.FileIndex}]))'*7500, 'b', 'LineWidth', 3);
end

handles.ASSLPickBoutsMainFig.ZoomSpecAxisLimits = handles.ASSLPickBoutsMainFig.SpecAxisLimits;
handles.ASSLPickBoutsMainFig.ZoomAmpAxisLimits = handles.ASSLPickBoutsMainFig.AmpAxisLimits;
handles.ASSLPickBoutsMainFig.ZoomLabelAxisLimits = handles.ASSLPickBoutsMainFig.LabelAxisLimits;

setappdata(findobj('Tag', 'ASSLPickBouts'), 'Data', handles.ASSLPickBoutsMainFig);
guidata(hObject, handles);

UpdateBoutNumber(handles);


% --- Executes on button press in DelAllBoutEdgesButton.
function DelAllBoutEdgesButton_Callback(hObject, eventdata, handles)
% hObject    handle to DelAllBoutEdgesButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

for i = 1:length(handles.ASSLPickBoutsMainFig.FileName),
    handles.BoutOnsetsOffsets{i} = [];
end
disp('Deleted all bout edges');
guidata(hObject, handles);



function InterBoutIntervalEdit_Callback(hObject, eventdata, handles)
% hObject    handle to InterBoutIntervalEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of InterBoutIntervalEdit as text
%        str2double(get(hObject,'String')) returns contents of InterBoutIntervalEdit as a double
handles.ASSLPickBoutsMainFig.InterBoutInterval = str2double(get(handles.InterBoutIntervalEdit, 'String'))*1000;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function InterBoutIntervalEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to InterBoutIntervalEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function BoutPaddingTimeEdit_Callback(hObject, eventdata, handles)
% hObject    handle to BoutPaddingTimeEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of BoutPaddingTimeEdit as text
%        str2double(get(hObject,'String')) returns contents of BoutPaddingTimeEdit as a double

handles.ASSLPickBoutsMainFig.BoutPaddingTime = str2double(get(hObject, 'String'))*1000;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function BoutPaddingTimeEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BoutPaddingTimeEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in LoResHiResToggle.
function LoResHiResToggle_Callback(hObject, eventdata, handles)
% hObject    handle to LoResHiResToggle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of LoResHiResToggle
handles.ASSLPickBoutsMainFig.LoResHiRes = get(hObject, 'Value');
if (handles.ASSLPickBoutsMainFig.LoResHiRes == 1)
    set(handles.LoResHiResToggle, 'String', 'High Res Spectrogram');
else
    set(handles.LoResHiResToggle, 'String', 'Low Res Spectrogram');
end
guidata(hObject, handles);

% --- Executes on button press in PlayFileButton.
function PlayFileButton_Callback(hObject, eventdata, handles)
% hObject    handle to PlayFileButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[RawData, Fs] = ASSLGetRawData(handles.ASSLPickBoutsMainFig.DirName, handles.ASSLPickBoutsMainFig.FileName{handles.ASSLPickBoutsMainFig.FileIndex}, handles.ASSLPickBoutsMainFig.FileType, handles.ASSLPickBoutsMainFig.SongChanNo);
RawData = RawData(ceil(handles.ASSLPickBoutsMainFig.ZoomSpecAxisLimits(1)*Fs):floor(handles.ASSLPickBoutsMainFig.ZoomSpecAxisLimits(2)*Fs));
soundsc(RawData, Fs);


% --- Executes on button press in ReCalcBoutEdgesButton.
function ReCalcBoutEdgesButton_Callback(hObject, eventdata, handles)
% hObject    handle to ReCalcBoutEdgesButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


for i = 1:length(handles.ASSLPickBoutsMainFig.FileName),
    [handles.BoutOnsetsOffsets{i}] = ASSLGetBoutEdges(handles.ASSLPickBoutsMainFig.SyllOnsets{i}, handles.ASSLPickBoutsMainFig.SyllOffsets{i}, handles.ASSLPickBoutsMainFig.InterBoutInterval, handles.ASSLPickBoutsMainFig.BoutPaddingTime, handles.ASSLPickBoutsMainFig.FileDur{i});
end

handles.ASSLPickBoutsMainFig.FileIndex = 1;
set(handles.SongFileNameTextLabel, 'String', ['Song File Name : ', handles.ASSLPickBoutsMainFig.FileName{handles.ASSLPickBoutsMainFig.FileIndex}, ' : #', num2str(handles.ASSLPickBoutsMainFig.FileIndex), ' of ', num2str(length(handles.ASSLPickBoutsMainFig.FileName)), ' files']);

[RawData, Fs] = ASSLGetRawData(handles.ASSLPickBoutsMainFig.DirName, handles.ASSLPickBoutsMainFig.FileName{handles.ASSLPickBoutsMainFig.FileIndex}, handles.ASSLPickBoutsMainFig.FileType, handles.ASSLPickBoutsMainFig.SongChanNo);

Time = (1:1:length(RawData))/Fs;

[LogAmplitude] = ASSLCalculateLogAmplitude(RawData, Fs, Time, handles.ASSLPickBoutsMainFig.FFTWinSizeSegmenting, handles.ASSLPickBoutsMainFig.FFTWinOverlapSegmenting);

%    set(handles.SongFileNameText, 'String', ['Song File Name : ', handles.ASSLPickBoutsMainFig.FileName{handles.ASSLPickBoutsMainFig.FileIndex}]);

[handles.ASSLPickBoutsMainFig.SpecAxisLimits, handles.ASSLPickBoutsMainFig.LabelAxisLimits, handles.ASSLPickBoutsMainFig.AmpAxisLimits] = ASSLReviewTMPlotData(handles, Time, LogAmplitude);

if (length(handles.ASSLPickBoutsMainFig.SyllOnsets{handles.ASSLPickBoutsMainFig.FileIndex}) > 1)
    axes(handles.ReviewSpecAxis);
    handles.ASSLPickBoutsMainFig.BoutPlotLine = plot([handles.BoutOnsetsOffsets{handles.ASSLPickBoutsMainFig.FileIndex}]'/1000, ones(size([handles.BoutOnsetsOffsets{handles.ASSLPickBoutsMainFig.FileIndex}]))'*7500, 'b', 'LineWidth', 3);
end

handles.ASSLPickBoutsMainFig.ZoomSpecAxisLimits = handles.ASSLPickBoutsMainFig.SpecAxisLimits;
handles.ASSLPickBoutsMainFig.ZoomAmpAxisLimits = handles.ASSLPickBoutsMainFig.AmpAxisLimits;
handles.ASSLPickBoutsMainFig.ZoomLabelAxisLimits = handles.ASSLPickBoutsMainFig.LabelAxisLimits;

% Update handles structure
guidata(hObject, handles);

UpdateBoutNumber(handles);

% Function to update the text string with total number of bouts
function [] = UpdateBoutNumber(handles)

TotalBoutNum = cell2mat(cellfun(@size, handles.BoutOnsetsOffsets, 'UniformOutput', 0));
TotalBoutNum = sum(TotalBoutNum(1:2:length(TotalBoutNum)));

TextString = get(handles.SongFileNameTextLabel, 'String');
TextString = [TextString, '; ', num2str(TotalBoutNum), ' bouts'];
set(handles.SongFileNameTextLabel, 'String', TextString);
