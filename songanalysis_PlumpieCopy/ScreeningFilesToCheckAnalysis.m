function varargout = ScreeningFilesToCheckAnalysis(varargin)
% SCREENINGFILESTOCHECKANALYSIS MATLAB code for ScreeningFilesToCheckAnalysis.fig
%      SCREENINGFILESTOCHECKANALYSIS, by itself, creates a new SCREENINGFILESTOCHECKANALYSIS or raises the existing
%      singleton*.
%
%      H = SCREENINGFILESTOCHECKANALYSIS returns the handle to a new SCREENINGFILESTOCHECKANALYSIS or the handle to
%      the existing singleton*.
%
%      SCREENINGFILESTOCHECKANALYSIS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SCREENINGFILESTOCHECKANALYSIS.M with the given input arguments.
%
%      SCREENINGFILESTOCHECKANALYSIS('Property','Value',...) creates a new SCREENINGFILESTOCHECKANALYSIS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ScreeningFilesToCheckAnalysis_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ScreeningFilesToCheckAnalysis_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ScreeningFilesToCheckAnalysis

% Last Modified by GUIDE v2.5 02-May-2019 02:31:25

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ScreeningFilesToCheckAnalysis_OpeningFcn, ...
                   'gui_OutputFcn',  @ScreeningFilesToCheckAnalysis_OutputFcn, ...
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


% --- Executes just before ScreeningFilesToCheckAnalysis is made visible.
function ScreeningFilesToCheckAnalysis_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ScreeningFilesToCheckAnalysis (see VARARGIN)

% Choose default command line output for ScreeningFilesToCheckAnalysis
handles.output = hObject;

% Some variables that need to be defined
handles.SFCA.ZoomStep = 6; % in seconds
set(handles.ZoomStepEdit, 'String', num2str(handles.SFCA.ZoomStep));
set(handles.Next5SecButton, 'String', ['Next ', num2str(handles.SFCA.ZoomStep), ' sec']);
set(handles.Prev5SecButton, 'String', ['Prev ', num2str(handles.SFCA.ZoomStep), ' sec']);

handles.SFCA.InterBoutInterval = 2; % in seconds
set(handles.InterBoutIntervalEdit, 'String', num2str(handles.SFCA.InterBoutInterval));

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ScreeningFilesToCheckAnalysis wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ScreeningFilesToCheckAnalysis_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in DataDirButton.
function DataDirButton_Callback(hObject, eventdata, handles)
% hObject    handle to DataDirButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.SFCA.DirName = uigetdir(pwd, 'Choose raw data directory');
cd(handles.SFCA.DirName);
set(handles.DataDirText, 'String', handles.SFCA.DirName);
guidata(hObject, handles);

% --- Executes on button press in SongFileListButton.
function SongFileListButton_Callback(hObject, eventdata, handles)
% hObject    handle to SongFileListButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[handles.SFCA.SongFileListName, handles.SFCA.SongFileListDirName] = uigetfile('*txt', 'Choose song file list');
set(handles.SongFileListText, 'String', handles.SFCA.SongFileListName);

% Load up the song file names
Fid = fopen(fullfile(handles.SFCA.SongFileListDirName, handles.SFCA.SongFileListName), 'r');
handles.SFCA.SongFiles = textscan(Fid, '%s', 'DeLimiter', '\n');
handles.SFCA.SongFiles = handles.SFCA.SongFiles{1};
fclose(Fid);

% Load up the note files and also get file length for each file
for i = 1:length(handles.SFCA.SongFiles),
    [RawData, Fs] = GetData(handles.SFCA.DirName, handles.SFCA.SongFiles{i}, handles.SFCA.FileTypes{handles.SFCA.ChosenFileType}, 0);
    handles.SFCA.FileLen(i) = length(RawData)/Fs;
    handles.SFCA.NoteInfo{i} = load(fullfile(handles.SFCA.DirName, 'ASSLNoteFiles', [handles.SFCA.SongFiles{i}, '.not.mat']));
end

handles.SFCA.SongFilesIndex = 1;
handles = SFCAPlotData(handles);

guidata(hObject, handles);


% --- Executes on selection change in FileTypeMenu.
function FileTypeMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileTypeMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns FileTypeMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from FileTypeMenu
handles.SFCA.FileTypes = cellstr(get(hObject, 'String'));
handles.SFCA.ChosenFileType = get(hObject, 'Value');
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function FileTypeMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FileTypeMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in PickVideoLogFileButton.
function PickVideoLogFileButton_Callback(hObject, eventdata, handles)
% hObject    handle to PickVideoLogFileButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in NextFileButton.
function NextFileButton_Callback(hObject, eventdata, handles)
% hObject    handle to NextFileButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.SFCA.SongFilesIndex = handles.SFCA.SongFilesIndex + 1;
if (handles.SFCA.SongFilesIndex > length(handles.SFCA.SongFiles))
    handles.SFCA.SongFilesIndex = length(handles.SFCA.SongFiles);
end
handles = SFCAPlotData(handles);
guidata(hObject, handles);

% --- Executes on button press in PrevFileButton.
function PrevFileButton_Callback(hObject, eventdata, handles)
% hObject    handle to PrevFileButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.SFCA.SongFilesIndex = handles.SFCA.SongFilesIndex - 1;
if (handles.SFCA.SongFilesIndex < 1)
    handles.SFCA.SongFilesIndex = 1;
end
handles = SFCAPlotData(handles);
guidata(hObject, handles);

% --- Executes on button press in Prev5SecButton.
function Prev5SecButton_Callback(hObject, eventdata, handles)
% hObject    handle to Prev5SecButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (diff(handles.SFCA.PlotAxisLimits(1:2)) <= handles.SFCA.ZoomStep)
    handles.SFCA.PlotAxisLimits(1:2) = [(handles.SFCA.PlotAxisLimits(1)*1.05 - handles.SFCA.ZoomStep) (handles.SFCA.PlotAxisLimits(1)*1.05)]; 
else
    handles.SFCA.PlotAxisLimits(1:2) = [0 handles.SFCA.ZoomStep]; 
end

if (handles.SFCA.PlotAxisLimits(1) < 0)
    handles.SFCA.PlotAxisLimits(1:2) = [0 handles.SFCA.ZoomStep];
end

handles.SFCA.LabelAxisLimits(1:2) = handles.SFCA.PlotAxisLimits(1:2);

axes(handles.SFCA_Axes);
axis(handles.SFCA.PlotAxisLimits);
axes(handles.SFCA_LabelAxes);
axis(handles.SFCA.LabelAxisLimits);
guidata(hObject, handles);

% --- Executes on button press in Next5SecButton.
function Next5SecButton_Callback(hObject, eventdata, handles)
% hObject    handle to Next5SecButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if (diff(handles.SFCA.PlotAxisLimits(1:2)) <= handles.SFCA.ZoomStep)
    handles.SFCA.PlotAxisLimits(1:2) = [(handles.SFCA.PlotAxisLimits(2)*0.95) (handles.SFCA.ZoomStep + handles.SFCA.PlotAxisLimits(2)*0.95)]; 
else
    handles.SFCA.PlotAxisLimits(1:2) = [0 handles.SFCA.ZoomStep]; 
end

if (handles.SFCA.PlotAxisLimits(2) > handles.SFCA.FileLen(handles.SFCA.SongFilesIndex))
    handles.SFCA.PlotAxisLimits(1:2) = handles.SFCA.FileLen(handles.SFCA.SongFilesIndex) - [handles.SFCA.ZoomStep 0];
end

handles.SFCA.LabelAxisLimits(1:2) = handles.SFCA.PlotAxisLimits(1:2);

axes(handles.SFCA_Axes);
axis(handles.SFCA.PlotAxisLimits);
axes(handles.SFCA_LabelAxes);
axis(handles.SFCA.LabelAxisLimits);
guidata(hObject, handles);


function ZoomStepEdit_Callback(hObject, eventdata, handles)
% hObject    handle to ZoomStepEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ZoomStepEdit as text
%        str2double(get(hObject,'String')) returns contents of ZoomStepEdit as a double
handles.SFCA.ZoomStep = str2double(get(hObject, 'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function ZoomStepEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ZoomStepEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function InterBoutIntervalEdit_Callback(hObject, eventdata, handles)
% hObject    handle to InterBoutIntervalEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of InterBoutIntervalEdit as text
%        str2double(get(hObject,'String')) returns contents of InterBoutIntervalEdit as a double
handles.SFCA.InterBoutInterval = str2double(get(hObject, 'String'));
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



function InterINIntervalEdit_Callback(hObject, eventdata, handles)
% hObject    handle to InterINIntervalEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of InterINIntervalEdit as text
%        str2double(get(hObject,'String')) returns contents of InterINIntervalEdit as a double
handles.SFCA.InterINInterval = str2double(get(hObject, 'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function InterINIntervalEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to InterINIntervalEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function MotifLabelsEdit_Callback(hObject, eventdata, handles)
% hObject    handle to MotifLabelsEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MotifLabelsEdit as text
%        str2double(get(hObject,'String')) returns contents of MotifLabelsEdit as a double
handles.SFCA.MotifLabels = get(hObject, 'String');
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function MotifLabelsEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MotifLabelsEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function INLabelsEdit_Callback(hObject, eventdata, handles)
% hObject    handle to INLabelsEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of INLabelsEdit as text
%        str2double(get(hObject,'String')) returns contents of INLabelsEdit as a double
handles.SFCA.INLabels = get(hObject, 'String');
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function INLabelsEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to INLabelsEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
