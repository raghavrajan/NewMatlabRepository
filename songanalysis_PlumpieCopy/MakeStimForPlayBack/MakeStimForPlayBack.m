function varargout = MakeStimForPlayBack(varargin)
% MAKESTIMFORPLAYBACK MATLAB code for MakeStimForPlayBack.fig
%      MAKESTIMFORPLAYBACK, by itself, creates a new MAKESTIMFORPLAYBACK or raises the existing
%      singleton*.
%
%      H = MAKESTIMFORPLAYBACK returns the handle to a new MAKESTIMFORPLAYBACK or the handle to
%      the existing singleton*.
%
%      MAKESTIMFORPLAYBACK('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MAKESTIMFORPLAYBACK.M with the given input arguments.
%
%      MAKESTIMFORPLAYBACK('Property','Value',...) creates a new MAKESTIMFORPLAYBACK or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MakeStimForPlayBack_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MakeStimForPlayBack_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MakeStimForPlayBack

% Last Modified by GUIDE v2.5 07-Apr-2014 11:44:32

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MakeStimForPlayBack_OpeningFcn, ...
                   'gui_OutputFcn',  @MakeStimForPlayBack_OutputFcn, ...
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


% --- Executes just before MakeStimForPlayBack is made visible.
function MakeStimForPlayBack_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MakeStimForPlayBack (see VARARGIN)

% Choose default command line output for MakeStimForPlayBack
handles.output = hObject;

handles.MSPB.StimDur = 15;
set(handles.StimDurationEdit, 'String', num2str(handles.MSPB.StimDur));

handles.MSPB.StimChoices = [{'Songs playing simultaneously on both speakers'}; {'Songs playing on only one speaker at a time'}];
handles.MSPB.StimChoice = 1;
set(handles.StimChoiceMenu, 'String', handles.MSPB.StimChoices);
set(handles.StimChoiceMenu, 'Value', handles.MSPB.StimChoice);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MakeStimForPlayBack wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = MakeStimForPlayBack_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in LoadWavFileButton.
function LoadWavFileButton_Callback(hObject, eventdata, handles)
% hObject    handle to LoadWavFileButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName, PathName] = uigetfile('*.wav');
MSPB_TrimWriteSongs(PathName, FileName);

% --- Executes on button press in WriteStimFileButton.
function WriteStimFileButton_Callback(hObject, eventdata, handles)
% hObject    handle to WriteStimFileButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.MSPB.SilenceLength = round(handles.MSPB.MaxStimLength + 1);
SilenceLen = round(handles.MSPB.SilenceFs * handles.MSPB.SilenceLength);
handles.MSPB.ActualSilence = repmat(handles.MSPB.Silence, floor(SilenceLen/size(handles.MSPB.Silence(:),1)), 1);
handles.MSPB.ActualSilence = [handles.MSPB.ActualSilence; handles.MSPB.Silence(1:(SilenceLen - size(handles.MSPB.ActualSilence(:),1)))];

SilenceLen = length(handles.MSPB.EqualLenStimuli{1});
handles.MSPB.EqualLenSilence = repmat(handles.MSPB.Silence, floor(SilenceLen/size(handles.MSPB.Silence(:),1)), 1);
handles.MSPB.EqualLenSilence = [handles.MSPB.EqualLenSilence; handles.MSPB.Silence(1:(SilenceLen - size(handles.MSPB.EqualLenSilence(:),1)))];

SongStimulus = [];

if (handles.MSPB.StimChoice == 1)
    Fid = fopen([handles.MSPB.SongFileListDirName, handles.MSPB.SongFileListName, '.Simultaneous.StimList.log'], 'w');
else
    Fid = fopen([handles.MSPB.SongFileListDirName, handles.MSPB.SongFileListName, '.Sequential.StimList.log'], 'w');
end

PrevLeftOrRight = [0 0 0];
Index = 0;
while (size(SongStimulus, 1) < (handles.MSPB.StimDur * 60 * handles.MSPB.StimuliFs{1}))
    if (handles.MSPB.StimChoice == 1)
        LeftStimulus = ceil(rand * length(handles.MSPB.EqualLenStimuli));
        RightStimulus = ceil(rand * length(handles.MSPB.EqualLenStimuli));
        SongStimulus = [SongStimulus; [[handles.MSPB.ActualSilence handles.MSPB.ActualSilence]; [handles.MSPB.EqualLenStimuli{LeftStimulus} handles.MSPB.EqualLenStimuli{RightStimulus}]]];
        fprintf(Fid, 'Left - #%i\tRight - #%i\n', LeftStimulus, RightStimulus);
    else
        if (mean(PrevLeftOrRight) == 1)
            LeftOrRight = 2;
        else
            if (mean(PrevLeftOrRight) == 2)
                LeftOrRight = 1;
            else
                LeftOrRight = ceil(rand*2);
            end
        end
        PrevLeftOrRight(mod(Index, 3) + 1) = LeftOrRight;
        Stimulus = ceil(rand * length(handles.MSPB.EqualLenStimuli));
        if (LeftOrRight == 1)
            SongStimulus = [SongStimulus; [[handles.MSPB.ActualSilence handles.MSPB.ActualSilence]; [handles.MSPB.EqualLenStimuli{Stimulus} handles.MSPB.EqualLenSilence]]];
            fprintf(Fid, 'Left - #%i\tRight - none\n', Stimulus);
        else
            SongStimulus = [SongStimulus; [[handles.MSPB.ActualSilence handles.MSPB.ActualSilence]; [handles.MSPB.EqualLenSilence handles.MSPB.EqualLenStimuli{Stimulus}]]];
            fprintf(Fid, 'Left - none\tRight - #%i\n', Stimulus);
        end
    end
    Index = Index + 1;
end
fclose(Fid);
if (handles.MSPB.StimChoice == 1)
    OutputFile = [handles.MSPB.SongFileListDirName, handles.MSPB.SongFileListName, '.Simultaneous.Stim.wav'];
else
    OutputFile = [handles.MSPB.SongFileListDirName, handles.MSPB.SongFileListName, '.Sequential.Stim.wav'];
end

wavwrite(SongStimulus, handles.MSPB.StimuliFs{1}, 16, OutputFile);
disp(['Finished writing stimulus to ', OutputFile]);

% --- Executes on button press in LoadStimuliFileListButton.
function LoadStimuliFileListButton_Callback(hObject, eventdata, handles)
% hObject    handle to LoadStimuliFileListButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

FileSep = filesep;

[handles.MSPB.SongFileListName, handles.MSPB.SongFileListDirName] = uigetfile('*', 'Choose input file list');
   
Fid = fopen([handles.MSPB.SongFileListDirName, handles.MSPB.SongFileListName], 'r');
TempFiles = textscan(Fid, '%s', 'DeLimiter', '\n');
handles.MSPB.SongFileNames = TempFiles{1};
fclose(Fid);

for i = 1:length(handles.MSPB.SongFileNames),
    SlashIndex = find((handles.MSPB.SongFileNames{i} == FileSep));
    if (~isempty(SlashIndex))
        handles.MSPB.SongFileNames{i} = handles.MSPB.SongFileNames{i}(SlashIndex(end)+1:end);
    end

    [handles.MSPB.Stimuli{i}, handles.MSPB.StimuliFs{i}] = ASSLGetRawData(handles.MSPB.SongFileListDirName, handles.MSPB.SongFileNames{i}, 'wav', 1);
    handles.MSPB.StimLength(i) = length(handles.MSPB.Stimuli{i})/handles.MSPB.StimuliFs{i};
end
handles.MSPB.MaxStimLength = max(handles.MSPB.StimLength);

for i = 1:length(handles.MSPB.Stimuli),
    handles.MSPB.EqualLenStimuli{i} = MSPB_MakeEqualLengthLeftRightStimuli(handles.MSPB.Stimuli{i}, handles.MSPB.Silence, handles.MSPB.SilenceFs, handles.MSPB.MaxStimLength);
end
disp('Finished loading stimuli');

guidata(hObject, handles);

% --- Executes on button press in LoadSilenceFileButton.
function LoadSilenceFileButton_Callback(hObject, eventdata, handles)
% hObject    handle to LoadSilenceFileButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[handles.MSPB.SilenceFileName, handles.MSPB.SilenceDirName] = uigetfile('*.silence.wav');
[handles.MSPB.Silence, handles.MSPB.SilenceFs] = ASSLGetRawData(handles.MSPB.SilenceDirName, handles.MSPB.SilenceFileName, 'wav', 1);
guidata(hObject, handles);

function StimDurationEdit_Callback(hObject, eventdata, handles)
% hObject    handle to StimDurationEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of StimDurationEdit as text
%        str2double(get(hObject,'String')) returns contents of StimDurationEdit as a double

handles.MSPB.StimDur = str2double(get(hObject, 'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function StimDurationEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to StimDurationEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in StimChoiceMenu.
function StimChoiceMenu_Callback(hObject, eventdata, handles)
% hObject    handle to StimChoiceMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns StimChoiceMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from StimChoiceMenu
handles.MSPB.StimChoice = get(hObject, 'Value');

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function StimChoiceMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to StimChoiceMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
