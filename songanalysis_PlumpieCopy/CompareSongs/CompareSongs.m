function varargout = CompareSongs(varargin)
% COMPARESONGS MATLAB code for CompareSongs.fig
%      COMPARESONGS, by itself, creates a new COMPARESONGS or raises the existing
%      singleton*.
%
%      H = COMPARESONGS returns the handle to a new COMPARESONGS or the handle to
%      the existing singleton*.
%
%      COMPARESONGS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in COMPARESONGS.M with the given input arguments.
%
%      COMPARESONGS('Property','Value',...) creates a new COMPARESONGS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before CompareSongs_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to CompareSongs_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help CompareSongs

% Last Modified by GUIDE v2.5 02-May-2014 13:33:13

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @CompareSongs_OpeningFcn, ...
                   'gui_OutputFcn',  @CompareSongs_OutputFcn, ...
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


% --- Executes just before CompareSongs is made visible.
function CompareSongs_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to CompareSongs (see VARARGIN)

% Choose default command line output for CompareSongs
handles.output = hObject;

handles.CS.MotifSyllLabels = 'abcde';
set(handles.MotifSyllLabelsEdit, 'String', handles.CS.MotifSyllLabels);

handles.CS.INLabels = 'i';
set(handles.INLabelsEdit, 'String', handles.CS.INLabels);

handles.CS.MotifInitiationSyllLabels = 'a';
set(handles.MotifInitiationSyllLabelsEdit, 'String', handles.CS.MotifInitiationSyllLabels);

handles.CS.MotifLabelsColor = [1 0 0];
Temp(:,:,1) = ones(10,100)*handles.CS.MotifLabelsColor(1);
Temp(:,:,2) = ones(10,100)*handles.CS.MotifLabelsColor(2);
Temp(:,:,3) = ones(10,100)*handles.CS.MotifLabelsColor(3);
set(handles.MotifLabelColorButton, 'CData', Temp);

handles.CS.INLabelsColor = [0 0 1];
Temp(:,:,1) = ones(10,100)*handles.CS.INLabelsColor(1);
Temp(:,:,2) = ones(10,100)*handles.CS.INLabelsColor(2);
Temp(:,:,3) = ones(10,100)*handles.CS.INLabelsColor(3);
set(handles.INLabelColorButton, 'CData', Temp);

handles.CS.SilenceColor = [1 1 1];
Temp(:,:,1) = ones(10,100)*handles.CS.SilenceColor(1);
Temp(:,:,2) = ones(10,100)*handles.CS.SilenceColor(2);
Temp(:,:,3) = ones(10,100)*handles.CS.SilenceColor(3);
set(handles.SilenceColorButton, 'CData', Temp);

handles.CS.OtherSyllsColor = [0 0 0];
Temp(:,:,1) = ones(10,100)*handles.CS.OtherSyllsColor(1);
Temp(:,:,2) = ones(10,100)*handles.CS.OtherSyllsColor(2);
Temp(:,:,3) = ones(10,100)*handles.CS.OtherSyllsColor(3);
set(handles.OtherSyllColorButton, 'CData', Temp);

handles.CS.NoDataColor = [0.8 0.8 0.8];
Temp(:,:,1) = ones(10,100)*handles.CS.NoDataColor(1);
Temp(:,:,2) = ones(10,100)*handles.CS.NoDataColor(2);
Temp(:,:,3) = ones(10,100)*handles.CS.NoDataColor(3);
set(handles.NoDataColorButton, 'CData', Temp);

handles.CS.BoutDefinitions = [{'Each new file is a bout'}; {'Each new file is a bout + inter-bout interval (sec)'}; {'Continuous data + inter-bout interval (sec)'}];
handles.CS.BoutDefChoice = 1;
set(handles.BoutDefChoiceMenu, 'String', handles.CS.BoutDefinitions);
set(handles.BoutDefChoiceMenu, 'Value', handles.CS.BoutDefChoice);

handles.CS.InterBoutInterval = 2000;
set(handles.InterBoutIntervalEdit, 'String', num2str(handles.CS.InterBoutInterval));
set(handles.InterBoutIntervalEdit, 'Enable', 'off');

handles.CS.PresentDir = pwd;
FileSep = filesep;
ProgramDir = which('CompareSongs');
FileSepIndex = find(ProgramDir == FileSep);
handles.AnalysisProgramsDir = [ProgramDir(1:FileSepIndex(end)), 'CSAnalysisPrograms'];
cd(handles.AnalysisProgramsDir);

AnalysisFiles = dir('*.m');

for i = 1:length(AnalysisFiles),
    handles.CS.AnalysisFiles{i} = AnalysisFiles(i).name(1:end-2);
    ListBoxString{i} = handles.CS.AnalysisFiles{i};
end
set(handles.AnalysisScriptsListBox, 'String', ListBoxString);
cd(handles.CS.PresentDir);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes CompareSongs wait for user response (see UIRESUME)
% uiwait(handles.CompareSongsMainWindow);


% --- Outputs from this function are returned to the command line.
function varargout = CompareSongs_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in LoadDataButton.
function LoadDataButton_Callback(hObject, eventdata, handles)
% hObject    handle to LoadDataButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (isfield(handles.CS, 'DataFileName'))
    handles.CS = rmfield(handles.CS, 'DataFileName');
end

if (isfield(handles.CS, 'DirName'))
    handles.CS = rmfield(handles.CS, 'DirName');
end
if (isfield(handles.CS, 'Data'))
    handles.CS = rmfield(handles.CS, 'Data');
end
if (isfield(handles.CS, 'Bouts'))
    handles.CS = rmfield(handles.CS, 'Bouts');
end

FileSep = filesep;
NoofDays = inputdlg('Enter the number of days that need to be analyzed', 'Number of days');
handles.CS.NoofDays = str2double(NoofDays{1});

for i = 1:handles.CS.NoofDays,
    [handles.CS.DataFileName{i}, handles.CS.DirName{i}] = uigetfile('*.ASSLData.mat', ['Choose data file for Day #', num2str(i)]);
    if (handles.CS.DirName{i}(end) ~= FileSep)
        handles.CS.DirName{i}(end+1) = FileSep;
    end
    Temp = load([handles.CS.DirName{i}, handles.CS.DataFileName{i}]);
    handles.CS.Data{i} = Temp.handles.ASSL;
    [handles.CS.Bouts{i}, handles.CS.BoutLens{i}, handles.CS.AllOnsets{i}, handles.CS.AllOffsets{i}, handles.CS.AllLabels{i}, handles.CS.AllFeats{i}, handles.CS.AllRemovalIndices{i}] = CSIdentifyBouts(Temp.handles.ASSL, handles.CS.BoutDefinitions{handles.CS.BoutDefChoice}, handles.CS.InterBoutInterval);
end
disp(['Finished loading data from ', num2str(handles.CS.NoofDays), ' days']);
guidata(hObject, handles);


% --- Executes on selection change in AnalysisScriptsListBox.
function AnalysisScriptsListBox_Callback(hObject, eventdata, handles)
% hObject    handle to AnalysisScriptsListBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns AnalysisScriptsListBox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from AnalysisScriptsListBox
feval(handles.CS.AnalysisFiles{get(hObject, 'Value')}, handles.CS);

% --- Executes during object creation, after setting all properties.
function AnalysisScriptsListBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AnalysisScriptsListBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in MotifSyllLabelsEdit.
function MotifSyllLabelsEdit_Callback(hObject, eventdata, handles)
% hObject    handle to MotifSyllLabelsEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.CS.MotifSyllLabels = get(hObject, 'String');
guidata(hObject, handles);

% --- Executes on button press in INLabelsEdit.
function INLabelsEdit_Callback(hObject, eventdata, handles)
% hObject    handle to INLabelsEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.CS.INLabels = get(hObject, 'String');
guidata(hObject, handles);

% --- Executes on selection change in BoutDefChoiceMenu.
function BoutDefChoiceMenu_Callback(hObject, eventdata, handles)
% hObject    handle to BoutDefChoiceMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns BoutDefChoiceMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from BoutDefChoiceMenu

handles.CS.BoutDefChoice = get(hObject, 'Value');
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function BoutDefChoiceMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BoutDefChoiceMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
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
handles.CS.InterBoutInterval = str2double(get(hObject, 'String'));
guidata(hObject, handles);


% --- Executes on button press in MotifLabelColorButton.
function MotifLabelColorButton_Callback(hObject, eventdata, handles)
% hObject    handle to MotifLabelColorButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.CS.MotifLabelsColor = uisetcolor;
Temp(:,:,1) = ones(10,100)*handles.CS.MotifLabelsColor(1);
Temp(:,:,2) = ones(10,100)*handles.CS.MotifLabelsColor(2);
Temp(:,:,3) = ones(10,100)*handles.CS.MotifLabelsColor(3);
set(handles.MotifLabelColorButton, 'CData', Temp);

guidata(hObject, handles);


% --- Executes on button press in INLabelColorButton.
function INLabelColorButton_Callback(hObject, eventdata, handles)
% hObject    handle to INLabelColorButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.CS.INLabelsColor = uisetcolor;
Temp(:,:,1) = ones(10,100)*handles.CS.INLabelsColor(1);
Temp(:,:,2) = ones(10,100)*handles.CS.INLabelsColor(2);
Temp(:,:,3) = ones(10,100)*handles.CS.INLabelsColor(3);
set(handles.INLabelColorButton, 'CData', Temp);
guidata(hObject, handles);


% --- Executes on button press in SilenceColorButton.
function SilenceColorButton_Callback(hObject, eventdata, handles)
% hObject    handle to SilenceColorButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.CS.SilenceColor = uisetcolor;
Temp(:,:,1) = ones(10,100)*handles.CS.SilenceColor(1);
Temp(:,:,2) = ones(10,100)*handles.CS.SilenceColor(2);
Temp(:,:,3) = ones(10,100)*handles.CS.SilenceColor(3);
set(handles.SilenceColorButton, 'CData', Temp);
guidata(hObject, handles);


% --- Executes on button press in OtherSyllColorButton.
function OtherSyllColorButton_Callback(hObject, eventdata, handles)
% hObject    handle to OtherSyllColorButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.CS.OtherSyllsColor = uisetcolor;
Temp(:,:,1) = ones(10,100)*handles.CS.OtherSyllsColor(1);
Temp(:,:,2) = ones(10,100)*handles.CS.OtherSyllsColor(2);
Temp(:,:,3) = ones(10,100)*handles.CS.OtherSyllsColor(3);
set(handles.OtherSyllColorButton, 'CData', Temp);
guidata(hObject, handles);


% --- Executes on button press in NoDataColorButton.
function NoDataColorButton_Callback(hObject, eventdata, handles)
% hObject    handle to NoDataColorButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.CS.NoDataColor = uisetcolor;
Temp(:,:,1) = ones(10,100)*handles.CS.NoDataColor(1);
Temp(:,:,2) = ones(10,100)*handles.CS.NoDataColor(2);
Temp(:,:,3) = ones(10,100)*handles.CS.NoDataColor(3);
set(handles.NoDataColorButton, 'CData', Temp);
guidata(hObject, handles);



function MotifInitiationSyllLabelsEdit_Callback(hObject, eventdata, handles)
% hObject    handle to MotifInitiationSyllLabelsEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MotifInitiationSyllLabelsEdit as text
%        str2double(get(hObject,'String')) returns contents of MotifInitiationSyllLabelsEdit as a double

handles.CS.MotifInitiationSyllLabels = get(hObject, 'String');
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function MotifInitiationSyllLabelsEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MotifInitiationSyllLabelsEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in CombineASSLDataFilesButton.
function CombineASSLDataFilesButton_Callback(hObject, eventdata, handles)
% hObject    handle to CombineASSLDataFilesButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

FileSep = filesep;
NoofFiles = inputdlg('Enter the number of files that need to be combined', 'Number of files');
NoofFiles = str2double(NoofFiles{1});

OutputFileName = [];

for i = 1:NoofFiles,
    [DataFileName{i}, DirName{i}] = uigetfile('*.ASSLData.mat', ['Choose data file for Day #', num2str(i)]);
    if (DirName{i}(end) ~= FileSep)
        DirName{i}(end+1) = FileSep;
    end
    Temp = load([DirName{i}, DataFileName{i}]);
    if (i == 1)
        ASSL.SyllOnsets = [];
        ASSL.SyllOffsets = [];
        ASSL.SyllLabels = [];
        ASSL.FileDur = [];
        for j = 1:length(Temp.handles.ASSL.ToBeUsedFeatures),
            eval(['ASSL.', Temp.handles.ASSL.ToBeUsedFeatures{j}, ' = [];']);
        end
        ASSL.ToBeUsedFeatures = Temp.handles.ASSL.ToBeUsedFeatures;
    end
    ASSL.SyllOnsets = [ASSL.SyllOnsets Temp.handles.ASSL.SyllOnsets];
    ASSL.SyllOffsets = [ASSL.SyllOffsets Temp.handles.ASSL.SyllOffsets];
    ASSL.SyllLabels = [ASSL.SyllLabels Temp.handles.ASSL.SyllLabels];
    ASSL.FileDur = [ASSL.FileDur Temp.handles.ASSL.FileDur];
    
    for j = 1:length(Temp.handles.ASSL.ToBeUsedFeatures),
        eval(['ASSL.', Temp.handles.ASSL.ToBeUsedFeatures{j}, ' = [ASSL.', Temp.handles.ASSL.ToBeUsedFeatures{j}, ' Temp.handles.ASSL.', Temp.handles.ASSL.ToBeUsedFeatures{j}, '];']);
    end
    
    OutputFileName = [OutputFileName, DataFileName{i}(1:(strfind(DataFileName{i}, '.ASSLData.mat') - 1)), '.'];
end

OutputFileName = [OutputFileName, 'ASSLData.mat'];
clear Temp;
handles.ASSL = ASSL;

OutputDir = uigetdir(pwd, 'Choose the directory for the combined file');
if (OutputDir(end) ~= FileSep)
    OutputDir(end+1) = FileSep;
end

save([OutputDir, OutputFileName], 'handles');
disp(['Finished combining data from ', num2str(NoofFiles), ' files']);