function varargout = ASSLPickBouts(varargin)
%ASSLPICKBOUTS M-file for ASSLPickBouts.fig
%      ASSLPICKBOUTS, by itself, creates a new ASSLPICKBOUTS or raises the existing
%      singleton*.
%
%      H = ASSLPICKBOUTS returns the handle to a new ASSLPICKBOUTS or the handle to
%      the existing singleton*.
%
%      ASSLPICKBOUTS('Property','Value',...) creates a new ASSLPICKBOUTS using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to ASSLPickBouts_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      ASSLPICKBOUTS('CALLBACK') and ASSLPICKBOUTS('CALLBACK',hObject,...) call the
%      local function named CALLBACK in ASSLPICKBOUTS.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ASSLPickBouts

% Last Modified by GUIDE v2.5 15-Oct-2014 09:51:55

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


% --- Executes just before ASSLPickBouts is made visible.
function ASSLPickBouts_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

handles.ASSLPickBouts.TimeStep = str2double(get(handles.TimeStepEdit, 'String'));
handles.ASSLPickBouts.InterBoutInterval = str2double(get(handles.InterBoutIntervalEdit, 'String'))*1000;
handles.ASSLPickBouts.BoutPaddingTime = str2double(get(handles.BoutPaddingTimeEdit, 'String'))*1000;

handles.ASSLPickBouts.LoResHiRes = 1;
set(handles.LoResHiResToggle, 'String', 'High Res Spectrogram');
set(handles.LoResHiResToggle, 'Value', handles.ASSLPickBouts.LoResHiRes);

if (nargin >= 1)
    handles.ASSLPickBouts = varargin{1};
    handles.ASSLPickBouts.LoResHiRes = 1;
    handles.ASSLPickBouts.TimeStep = str2double(get(handles.TimeStepEdit, 'String'));
    handles.ASSLPickBouts.InterBoutInterval = str2double(get(handles.InterBoutIntervalEdit, 'String'))*1000;
    handles.ASSLPickBouts.BoutPaddingTime = str2double(get(handles.BoutPaddingTimeEdit, 'String'))*1000;
    
    for i = 1:length(handles.ASSLPickBouts.FileName),
        [handles.BoutOnsetsOffsets{i}] = ASSLGetBoutEdges(handles.ASSLPickBouts.SyllOnsets{i}, handles.ASSLPickBouts.SyllOffsets{i}, handles.ASSLPickBouts.InterBoutInterval, handles.ASSLPickBouts.BoutPaddingTime, handles.ASSLPickBouts.FileDur{i});
    end
    
    handles.ASSLPickBouts.FileIndex = 1;
    set(handles.SongFileNameTextLabel, 'String', ['Song File Name : ', handles.ASSLPickBouts.FileName{handles.ASSLPickBouts.FileIndex}, ' : #', num2str(handles.ASSLPickBouts.FileIndex), ' of ', num2str(length(handles.ASSLPickBouts.FileName)), ' files']);
    
    [RawData, Fs] = ASSLGetRawData(handles.ASSLPickBouts.DirName, handles.ASSLPickBouts.FileName{handles.ASSLPickBouts.FileIndex}, handles.ASSLPickBouts.FileType, handles.ASSLPickBouts.SongChanNo);

    Time = (1:1:length(RawData))/Fs;

    [LogAmplitude] = ASSLCalculateLogAmplitude(RawData, Fs, Time, handles.ASSLPickBouts.FFTWinSizeSegmenting, handles.ASSLPickBouts.FFTWinOverlapSegmenting);

%    set(handles.SongFileNameText, 'String', ['Song File Name : ', handles.ASSLPickBouts.FileName{handles.ASSLPickBouts.FileIndex}]);

    [handles.ASSLPickBouts.SpecAxisLimits, handles.ASSLPickBouts.LabelAxisLimits, handles.ASSLPickBouts.AmpAxisLimits] = ASSLReviewTMPlotData(handles, Time, LogAmplitude);
    
    if (length(handles.ASSLPickBouts.SyllOnsets{handles.ASSLPickBouts.FileIndex}) > 1)
        axes(handles.ReviewSpecAxis);
        handles.ASSLPickBouts.BoutPlotLine = plot([handles.BoutOnsetsOffsets{handles.ASSLPickBouts.FileIndex}]'/1000, ones(size([handles.BoutOnsetsOffsets{handles.ASSLPickBouts.FileIndex}]))'*7500, 'b', 'LineWidth', 3);
    end
end

handles.ASSLPickBouts.ZoomSpecAxisLimits = handles.ASSLPickBouts.SpecAxisLimits;
handles.ASSLPickBouts.ZoomAmpAxisLimits = handles.ASSLPickBouts.AmpAxisLimits;
handles.ASSLPickBouts.ZoomLabelAxisLimits = handles.ASSLPickBouts.LabelAxisLimits;

% Update handles structure
guidata(hObject, handles);

UpdateBoutNumber(handles);

% UIWAIT makes ASSLPickBouts wait for user response (see UIRESUME)
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
for i = 1:length(handles.ASSLPickBouts.FileName),
    if (~isempty(handles.BoutOnsetsOffsets{i}))
        onsets = handles.BoutOnsetsOffsets{i}(:,1);
        offsets = handles.BoutOnsetsOffsets{i}(:,2);
        labels = repmat('B', 1, length(onsets));
    else
        onsets = [];
        offsets = [];
        labels = [];
    end
    min_int = handles.ASSLPickBouts.InterBoutInterval;
    save([handles.ASSLPickBouts.FileName{i}, '.not.mat'], 'onsets', 'offsets', 'labels', 'min_int');
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
handles.ASSLPickBouts.TimeStep = str2double(get(hObject, 'String'));
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
            TempBouts = handles.BoutOnsetsOffsets{handles.ASSLPickBouts.FileIndex};
            if (~isempty(TempBouts)),
                for i = 1:size(TempBouts,1),
                    if ((x(1) >= TempBouts(i,1)) && (x(1) <= TempBouts(i,2)))
                        TempBouts(i,:) = [];
                        break;
                    end
                end
                handles.BoutOnsetsOffsets{handles.ASSLPickBouts.FileIndex} = TempBouts;
                axes(handles.ReviewSpecAxis);
                delete(handles.ASSLPickBouts.BoutPlotLine);
                handles.ASSLPickBouts.BoutPlotLine = plot([handles.BoutOnsetsOffsets{handles.ASSLPickBouts.FileIndex}]'/1000, ones(size([handles.BoutOnsetsOffsets{handles.ASSLPickBouts.FileIndex}]))'*7500, 'b', 'LineWidth', 3);
                UpdateBoutNumber(handles);
            end
        
        case 78 % typing N
            NextTimeButton_Callback(hObject, eventdata, handles);
            handles.ASSLPickBouts = getappdata(findobj('Tag', 'ASSLPickBouts'), 'Data');

        case 80 % typing P
            PrevTimeButton_Callback(hObject, eventdata, handles);
            handles.ASSLPickBouts = getappdata(findobj('Tag', 'ASSLPickBouts'), 'Data');
            
        case 110 % typing n
            NextFileButton_Callback(hObject, eventdata, handles);
            handles.ASSLPickBouts = getappdata(findobj('Tag', 'ASSLPickBouts'), 'Data');
            
        case 112 % typing p
            PrevFileButton_Callback(hObject, eventdata, handles);
            handles.ASSLPickBouts = getappdata(findobj('Tag', 'ASSLPickBouts'), 'Data');
            
    end
end
guidata(hObject, handles);


% --- Executes on button press in NextFileButton.
function NextFileButton_Callback(hObject, eventdata, handles)
% hObject    handle to NextFileButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.ASSLPickBouts.FileIndex = handles.ASSLPickBouts.FileIndex + 1;
if (handles.ASSLPickBouts.FileIndex < 1)
    handles.ASSLPickBouts.FileIndex = 1;
end

if (handles.ASSLPickBouts.FileIndex > length(handles.ASSLPickBouts.FileName))
    handles.ASSLPickBouts.FileIndex = length(handles.ASSLPickBouts.FileName);
end

set(handles.SongFileNameTextLabel, 'String', ['Song File Name : ', handles.ASSLPickBouts.FileName{handles.ASSLPickBouts.FileIndex}, ' : #', num2str(handles.ASSLPickBouts.FileIndex), ' of ', num2str(length(handles.ASSLPickBouts.FileName)), ' files']);
[RawData, Fs] = ASSLGetRawData(handles.ASSLPickBouts.DirName, handles.ASSLPickBouts.FileName{handles.ASSLPickBouts.FileIndex}, handles.ASSLPickBouts.FileType, handles.ASSLPickBouts.SongChanNo);
Time = (1:1:length(RawData))/Fs;
[LogAmplitude] = ASSLCalculateLogAmplitude(RawData, Fs, Time, handles.ASSLPickBouts.FFTWinSizeSegmenting, handles.ASSLPickBouts.FFTWinOverlapSegmenting);

[handles.ASSLPickBouts.SpecAxisLimits, handles.ASSLPickBouts.LabelAxisLimits, handles.ASSLPickBouts.AmpAxisLimits] = ASSLReviewTMPlotData(handles, Time, LogAmplitude);

if (length(handles.ASSLPickBouts.SyllOnsets{handles.ASSLPickBouts.FileIndex}) > 1)
    axes(handles.ReviewSpecAxis);
    handles.ASSLPickBouts.BoutPlotLine = plot([handles.BoutOnsetsOffsets{handles.ASSLPickBouts.FileIndex}]'/1000, ones(size([handles.BoutOnsetsOffsets{handles.ASSLPickBouts.FileIndex}]))'*7500, 'b', 'LineWidth', 3);
end

handles.ASSLPickBouts.ZoomSpecAxisLimits = handles.ASSLPickBouts.SpecAxisLimits;
handles.ASSLPickBouts.ZoomAmpAxisLimits = handles.ASSLPickBouts.AmpAxisLimits;
handles.ASSLPickBouts.ZoomLabelAxisLimits = handles.ASSLPickBouts.LabelAxisLimits;

setappdata(findobj('Tag', 'ASSLPickBouts'), 'Data', handles.ASSLPickBouts);
guidata(hObject, handles);

UpdateBoutNumber(handles);

% --- Executes on button press in PrevFileButton.
function PrevFileButton_Callback(hObject, eventdata, handles)
% hObject    handle to PrevFileButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.ASSLPickBouts.FileIndex = handles.ASSLPickBouts.FileIndex - 1;
if (handles.ASSLPickBouts.FileIndex < 1)
    handles.ASSLPickBouts.FileIndex = 1;
end

if (handles.ASSLPickBouts.FileIndex > length(handles.ASSLPickBouts.FileName))
    handles.ASSLPickBouts.FileIndex = length(handles.ASSLPickBouts.FileName);
end

set(handles.SongFileNameTextLabel, 'String', ['Song File Name : ', handles.ASSLPickBouts.FileName{handles.ASSLPickBouts.FileIndex}, ' : #', num2str(handles.ASSLPickBouts.FileIndex), ' of ', num2str(length(handles.ASSLPickBouts.FileName)), ' files']);
[RawData, Fs] = ASSLGetRawData(handles.ASSLPickBouts.DirName, handles.ASSLPickBouts.FileName{handles.ASSLPickBouts.FileIndex}, handles.ASSLPickBouts.FileType, handles.ASSLPickBouts.SongChanNo);
Time = (1:1:length(RawData))/Fs;
[LogAmplitude] = ASSLCalculateLogAmplitude(RawData, Fs, Time, handles.ASSLPickBouts.FFTWinSizeSegmenting, handles.ASSLPickBouts.FFTWinOverlapSegmenting);

[handles.ASSLPickBouts.SpecAxisLimits, handles.ASSLPickBouts.LabelAxisLimits, handles.ASSLPickBouts.AmpAxisLimits] = ASSLReviewTMPlotData(handles, Time, LogAmplitude);

if (length(handles.ASSLPickBouts.SyllOnsets{handles.ASSLPickBouts.FileIndex}) > 1)
    axes(handles.ReviewSpecAxis);
    handles.ASSLPickBouts.BoutPlotLine = plot([handles.BoutOnsetsOffsets{handles.ASSLPickBouts.FileIndex}]'/1000, ones(size([handles.BoutOnsetsOffsets{handles.ASSLPickBouts.FileIndex}]))'*7500, 'b', 'LineWidth', 3);
end

handles.ASSLPickBouts.ZoomSpecAxisLimits = handles.ASSLPickBouts.SpecAxisLimits;
handles.ASSLPickBouts.ZoomAmpAxisLimits = handles.ASSLPickBouts.AmpAxisLimits;
handles.ASSLPickBouts.ZoomLabelAxisLimits = handles.ASSLPickBouts.LabelAxisLimits;

setappdata(findobj('Tag', 'ASSLPickBouts'), 'Data', handles.ASSLPickBouts);
guidata(hObject, handles);

UpdateBoutNumber(handles);

% --- Executes on button press in NextTimeButton.
function NextTimeButton_Callback(hObject, eventdata, handles)
% hObject    handle to NextTimeButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ((handles.ASSLPickBouts.ZoomSpecAxisLimits(2) - handles.ASSLPickBouts.ZoomSpecAxisLimits(1)) > 1.1*handles.ASSLPickBouts.TimeStep)
    handles.ASSLPickBouts.ZoomSpecAxisLimits(2) = handles.ASSLPickBouts.SpecAxisLimits(1) + handles.ASSLPickBouts.TimeStep;
else
    handles.ASSLPickBouts.ZoomSpecAxisLimits(1) = handles.ASSLPickBouts.ZoomSpecAxisLimits(1) + 0.9*handles.ASSLPickBouts.TimeStep;
    handles.ASSLPickBouts.ZoomSpecAxisLimits(2) = handles.ASSLPickBouts.ZoomSpecAxisLimits(1) + handles.ASSLPickBouts.TimeStep;
end

if (handles.ASSLPickBouts.ZoomSpecAxisLimits(1) < 0)
    handles.ASSLPickBouts.ZoomSpecAxisLimits(1) = 0;
    handles.ASSLPickBouts.ZoomSpecAxisLimits(2) = handles.ASSLPickBouts.ZoomSpecAxisLimits(1) + handles.ASSLPickBouts.TimeStep;
end

if (handles.ASSLPickBouts.ZoomSpecAxisLimits(2) > handles.ASSLPickBouts.SpecAxisLimits(2))
    handles.ASSLPickBouts.ZoomSpecAxisLimits(2) = handles.ASSLPickBouts.SpecAxisLimits(2) + handles.ASSLPickBouts.TimeStep*0.1;
    handles.ASSLPickBouts.ZoomSpecAxisLimits(1) = handles.ASSLPickBouts.ZoomSpecAxisLimits(2) - handles.ASSLPickBouts.TimeStep;
end

handles.ASSLPickBouts.ZoomAmpAxisLimits(1:2) = handles.ASSLPickBouts.ZoomSpecAxisLimits(1:2);
handles.ASSLPickBouts.ZoomLabelAxisLimits(1:2) = handles.ASSLPickBouts.ZoomSpecAxisLimits(1:2);

axes(handles.ReviewSpecAxis);
axis(handles.ASSLPickBouts.ZoomSpecAxisLimits);

axes(handles.ReviewLabelAxis);
axis(handles.ASSLPickBouts.ZoomLabelAxisLimits);

axes(handles.ReviewAmplitudeAxis);
axis(handles.ASSLPickBouts.ZoomAmpAxisLimits);

setappdata(findobj('Tag', 'ASSLPickBouts'), 'Data', handles.ASSLPickBouts);
guidata(hObject, handles);

% --- Executes on button press in PrevTimeButton.
function PrevTimeButton_Callback(hObject, eventdata, handles)
% hObject    handle to PrevTimeButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


if ((handles.ASSLPickBouts.ZoomSpecAxisLimits(2) - handles.ASSLPickBouts.ZoomSpecAxisLimits(1)) > 1.1*handles.ASSLPickBouts.TimeStep)
    handles.ASSLPickBouts.ZoomSpecAxisLimits(2) = handles.ASSLPickBouts.SpecAxisLimits(1) + handles.ASSLPickBouts.TimeStep;
else
    handles.ASSLPickBouts.ZoomSpecAxisLimits(1) = handles.ASSLPickBouts.ZoomSpecAxisLimits(1) - 0.9*handles.ASSLPickBouts.TimeStep;
    handles.ASSLPickBouts.ZoomSpecAxisLimits(2) = handles.ASSLPickBouts.ZoomSpecAxisLimits(1) + handles.ASSLPickBouts.TimeStep;
end

if (handles.ASSLPickBouts.ZoomSpecAxisLimits(1) < 0)
    handles.ASSLPickBouts.ZoomSpecAxisLimits(1) = 0;
    handles.ASSLPickBouts.ZoomSpecAxisLimits(2) = handles.ASSLPickBouts.ZoomSpecAxisLimits(1) + handles.ASSLPickBouts.TimeStep;
end

if (handles.ASSLPickBouts.ZoomSpecAxisLimits(2) > handles.ASSLPickBouts.SpecAxisLimits(2))
    handles.ASSLPickBouts.ZoomSpecAxisLimits(2) = handles.ASSLPickBouts.SpecAxisLimits(2) + handles.ASSLPickBouts.TimeStep*0.1;
    handles.ASSLPickBouts.ZoomSpecAxisLimits(1) = handles.ASSLPickBouts.ZoomSpecAxisLimits(2) - handles.ASSLPickBouts.TimeStep;
end


handles.ASSLPickBouts.ZoomAmpAxisLimits(1:2) = handles.ASSLPickBouts.ZoomSpecAxisLimits(1:2);
handles.ASSLPickBouts.ZoomLabelAxisLimits(1:2) = handles.ASSLPickBouts.ZoomSpecAxisLimits(1:2);

axes(handles.ReviewSpecAxis);
axis(handles.ASSLPickBouts.ZoomSpecAxisLimits);

axes(handles.ReviewLabelAxis);
axis(handles.ASSLPickBouts.ZoomLabelAxisLimits);

axes(handles.ReviewAmplitudeAxis);
axis(handles.ASSLPickBouts.ZoomAmpAxisLimits);

setappdata(findobj('Tag', 'ASSLPickBouts'), 'Data', handles.ASSLPickBouts);
guidata(hObject, handles);


% --- Executes on button press in JumpFileButton.
function JumpFileButton_Callback(hObject, eventdata, handles)
% hObject    handle to JumpFileButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

NewFileIndex = inputdlg('Enter the file # to jump to', 'File #');
NewFileIndex = str2double(NewFileIndex{1});

handles.ASSLPickBouts.FileIndex = NewFileIndex;
if (handles.ASSLPickBouts.FileIndex < 1)
    handles.ASSLPickBouts.FileIndex = 1;
end

if (handles.ASSLPickBouts.FileIndex > length(handles.ASSLPickBouts.FileName))
    handles.ASSLPickBouts.FileIndex = length(handles.ASSLPickBouts.FileName);
end

set(handles.SongFileNameTextLabel, 'String', ['Song File Name : ', handles.ASSLPickBouts.FileName{handles.ASSLPickBouts.FileIndex}, ' : #', num2str(handles.ASSLPickBouts.FileIndex), ' of ', num2str(length(handles.ASSLPickBouts.FileName)), ' files']);
[RawData, Fs] = ASSLGetRawData(handles.ASSLPickBouts.DirName, handles.ASSLPickBouts.FileName{handles.ASSLPickBouts.FileIndex}, handles.ASSLPickBouts.FileType, handles.ASSLPickBouts.SongChanNo);
Time = (1:1:length(RawData))/Fs;
[LogAmplitude] = ASSLCalculateLogAmplitude(RawData, Fs, Time, handles.ASSLPickBouts.FFTWinSizeSegmenting, handles.ASSLPickBouts.FFTWinOverlapSegmenting);

[handles.ASSLPickBouts.SpecAxisLimits, handles.ASSLPickBouts.LabelAxisLimits, handles.ASSLPickBouts.AmpAxisLimits] = ASSLReviewTMPlotData(handles, Time, LogAmplitude);

if (length(handles.ASSLPickBouts.SyllOnsets{handles.ASSLPickBouts.FileIndex}) > 1)
    axes(handles.ReviewSpecAxis);
    handles.ASSLPickBouts.BoutPlotLine = plot([handles.BoutOnsetsOffsets{handles.ASSLPickBouts.FileIndex}]'/1000, ones(size([handles.BoutOnsetsOffsets{handles.ASSLPickBouts.FileIndex}]))'*7500, 'b', 'LineWidth', 3);
end

handles.ASSLPickBouts.ZoomSpecAxisLimits = handles.ASSLPickBouts.SpecAxisLimits;
handles.ASSLPickBouts.ZoomAmpAxisLimits = handles.ASSLPickBouts.AmpAxisLimits;
handles.ASSLPickBouts.ZoomLabelAxisLimits = handles.ASSLPickBouts.LabelAxisLimits;

setappdata(findobj('Tag', 'ASSLPickBouts'), 'Data', handles.ASSLPickBouts);
guidata(hObject, handles);

UpdateBoutNumber(handles);


% --- Executes on button press in DelAllBoutEdgesButton.
function DelAllBoutEdgesButton_Callback(hObject, eventdata, handles)
% hObject    handle to DelAllBoutEdgesButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

for i = 1:length(handles.ASSLPickBouts.FileName),
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
handles.ASSLPickBouts.InterBoutInterval = str2double(get(handles.InterBoutIntervalEdit, 'String'))*1000;
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

handles.ASSLPickBouts.BoutPaddingTime = str2double(get(hObject, 'String'))*1000;
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
handles.ASSLPickBouts.LoResHiRes = get(hObject, 'Value');
if (handles.ASSLPickBouts.LoResHiRes == 1)
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

[RawData, Fs] = ASSLGetRawData(handles.ASSLPickBouts.DirName, handles.ASSLPickBouts.FileName{handles.ASSLPickBouts.FileIndex}, handles.ASSLPickBouts.FileType, handles.ASSLPickBouts.SongChanNo);
RawData = RawData(ceil(handles.ASSLPickBouts.ZoomSpecAxisLimits(1)*Fs):floor(handles.ASSLPickBouts.ZoomSpecAxisLimits(2)*Fs));
soundsc(RawData, Fs);


% --- Executes on button press in ReCalcBoutEdgesButton.
function ReCalcBoutEdgesButton_Callback(hObject, eventdata, handles)
% hObject    handle to ReCalcBoutEdgesButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


for i = 1:length(handles.ASSLPickBouts.FileName),
    [handles.BoutOnsetsOffsets{i}] = ASSLGetBoutEdges(handles.ASSLPickBouts.SyllOnsets{i}, handles.ASSLPickBouts.SyllOffsets{i}, handles.ASSLPickBouts.InterBoutInterval, handles.ASSLPickBouts.BoutPaddingTime, handles.ASSLPickBouts.FileDur{i});
end

handles.ASSLPickBouts.FileIndex = 1;
set(handles.SongFileNameTextLabel, 'String', ['Song File Name : ', handles.ASSLPickBouts.FileName{handles.ASSLPickBouts.FileIndex}, ' : #', num2str(handles.ASSLPickBouts.FileIndex), ' of ', num2str(length(handles.ASSLPickBouts.FileName)), ' files']);

[RawData, Fs] = ASSLGetRawData(handles.ASSLPickBouts.DirName, handles.ASSLPickBouts.FileName{handles.ASSLPickBouts.FileIndex}, handles.ASSLPickBouts.FileType, handles.ASSLPickBouts.SongChanNo);

Time = (1:1:length(RawData))/Fs;

[LogAmplitude] = ASSLCalculateLogAmplitude(RawData, Fs, Time, handles.ASSLPickBouts.FFTWinSizeSegmenting, handles.ASSLPickBouts.FFTWinOverlapSegmenting);

%    set(handles.SongFileNameText, 'String', ['Song File Name : ', handles.ASSLPickBouts.FileName{handles.ASSLPickBouts.FileIndex}]);

[handles.ASSLPickBouts.SpecAxisLimits, handles.ASSLPickBouts.LabelAxisLimits, handles.ASSLPickBouts.AmpAxisLimits] = ASSLReviewTMPlotData(handles, Time, LogAmplitude);

if (length(handles.ASSLPickBouts.SyllOnsets{handles.ASSLPickBouts.FileIndex}) > 1)
    axes(handles.ReviewSpecAxis);
    handles.ASSLPickBouts.BoutPlotLine = plot([handles.BoutOnsetsOffsets{handles.ASSLPickBouts.FileIndex}]'/1000, ones(size([handles.BoutOnsetsOffsets{handles.ASSLPickBouts.FileIndex}]))'*7500, 'b', 'LineWidth', 3);
end

handles.ASSLPickBouts.ZoomSpecAxisLimits = handles.ASSLPickBouts.SpecAxisLimits;
handles.ASSLPickBouts.ZoomAmpAxisLimits = handles.ASSLPickBouts.AmpAxisLimits;
handles.ASSLPickBouts.ZoomLabelAxisLimits = handles.ASSLPickBouts.LabelAxisLimits;

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
