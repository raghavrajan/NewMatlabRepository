function varargout = MicroLesionMainAnalysis(varargin)
% MICROLESIONMAINANALYSIS MATLAB code for MicroLesionMainAnalysis.fig
%      MICROLESIONMAINANALYSIS, by itself, creates a new MICROLESIONMAINANALYSIS or raises the existing
%      singleton*.
%
%      H = MICROLESIONMAINANALYSIS returns the handle to a new MICROLESIONMAINANALYSIS or the handle to
%      the existing singleton*.
%
%      MICROLESIONMAINANALYSIS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MICROLESIONMAINANALYSIS.M with the given input arguments.
%
%      MICROLESIONMAINANALYSIS('Property','Value',...) creates a new MICROLESIONMAINANALYSIS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MicroLesionMainAnalysis_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MicroLesionMainAnalysis_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MicroLesionMainAnalysis

% Last Modified by GUIDE v2.5 28-Nov-2013 16:19:49

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MicroLesionMainAnalysis_OpeningFcn, ...
                   'gui_OutputFcn',  @MicroLesionMainAnalysis_OutputFcn, ...
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


% --- Executes just before MicroLesionMainAnalysis is made visible.
function MicroLesionMainAnalysis_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MicroLesionMainAnalysis (see VARARGIN)

% Choose default command line output for MicroLesionMainAnalysis
handles.output = hObject;

% Initialize variable and gui data
handles.AllFileTypes = cellstr(get(handles.FileTypeMenu, 'String'));

if (~ispc)
    [status, result] = unix('whoami');
    handles.OutputDir = ['/home/', result(1:end-1), '/MicrolesionAnalysisResults/'];
    if (~exist(['/home/', result(1:end-1), '/MicrolesionAnalysisResults/'], 'dir'))
        eval(['!mkdir ''/home/' result(1:end-1) '/MicrolesionAnalysisResults''']);
    end
    set(handles.OutputDirEdit, 'String', handles.OutputDir);
end

PresentDir = pwd;

% Initialize axis to be blank without any axis
axes(handles.MLMATemplateAxis);
set(gca, 'XTick', []);
set(gca, 'YTick', []);
set(gca, 'Visible', 'off');

axes(handles.MLMASyllTemplateAxis);
set(gca, 'XTick', []);
set(gca, 'YTick', []);
set(gca, 'Visible', 'off');

axes(handles.ProgressAxis);
set(gca, 'XTick', []);
set(gca, 'YTick', []);
set(gca, 'Visible', 'off');

handles.SyllTemplates.SyllableTemplates = [];

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MicroLesionMainAnalysis wait for user response (see UIRESUME)
% uiwait(handles.MicroLesionAnalysisMainWindow);


% --- Outputs from this function are returned to the command line.
function varargout = MicroLesionMainAnalysis_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function BirdNameEdit_Callback(hObject, eventdata, handles)
% hObject    handle to BirdNameEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of BirdNameEdit as text
%        str2double(get(hObject,'String')) returns contents of BirdNameEdit as a double

handles.BirdName = get(hObject, 'String');
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function BirdNameEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BirdNameEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ChooseMotifTemplateFileButton.
function ChooseMotifTemplateFileButton_Callback(hObject, eventdata, handles)
% hObject    handle to ChooseMotifTemplateFileButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[TemplateFileName, TemplateFileDir] = uigetfile('*.mat', 'Choose a template file');
handles.MotifTemplate = load(fullfile(TemplateFileDir, TemplateFileName));

ZeroStretchIndex = find(([handles.MotifTemplate.MotifTemplate.TimeStretch] == 0) & ([handles.MotifTemplate.MotifTemplate.FreqStretch] == 0));

if (~isempty(ZeroStretchIndex))
    axes(handles.MLMATemplateAxis);
    contourf(handles.MotifTemplate.MotifTemplate(ZeroStretchIndex).MotifTemplate);
    colorbar;
    set(gca, 'YTick', []);
    set(gca, 'XTick', []);
    set(handles.MotifTemplateFileNameText, 'String', TemplateFileName);
end

guidata(hObject, handles);


% --- Executes on button press in ChooseSyllTemplatesButton.
function ChooseSyllTemplatesButton_Callback(hObject, eventdata, handles)
% hObject    handle to ChooseSyllTemplatesButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[TemplateFileName, TemplateFileDir] = uigetfile('*.mat', 'Choose a template file');
handles.SyllTemplates = load(fullfile(TemplateFileDir, TemplateFileName));

cla(handles.MLMASyllTemplateAxis);
axes(handles.MLMASyllTemplateAxis);
hold on;

Increment = 0;
for i = 1:length(handles.SyllTemplates.SyllableTemplates),
    ZeroStretchIndex = find(([handles.SyllTemplates.SyllableTemplates{i}{1}.MotifTemplate.TimeStretch] == 0) & ([handles.SyllTemplates.SyllableTemplates{i}{1}.MotifTemplate.FreqStretch] == 0));
    if (~isempty(ZeroStretchIndex))
        SyllTemplate = handles.SyllTemplates.SyllableTemplates{i}{1}.MotifTemplate(ZeroStretchIndex).MotifTemplate;
        contourf((1:1:size(SyllTemplate,2))+Increment, 1:1:size(SyllTemplate,1), SyllTemplate);
        
        SyllLabel = handles.SyllTemplates.SyllableTemplates{i}{1}.MotifTemplate(ZeroStretchIndex).Label;
        text(size(SyllTemplate,2)/2 + Increment, size(SyllTemplate,1) + 5, SyllLabel, 'FontWeight', 'bold', 'FontSize', 14);
                
        Increment = Increment + size(SyllTemplate,2) + 5;
        
        axis tight;
        set(gca, 'YTick', []);
        set(gca, 'XTick', []);
    end
end
colorbar;
axis tight;
temp = axis;
temp(4) = temp(4) + 10;
axis(temp);

set(handles.SyllTemplateFileNameText, 'String', TemplateFileName);
guidata(hObject, handles);

% --- Executes on selection change in FileTypeMenu.
function FileTypeMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileTypeMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns FileTypeMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from FileTypeMenu

handles.AllFileTypes = cellstr(get(hObject, 'String'));

if (strfind(handles.AllFileTypes{get(hObject, 'Value')}, 'Observer'))
    handles.FileType = 'obs';
else
    if (strfind(handles.AllFileTypes{get(hObject, 'Value')}, 'OKrank'))
        handles.FileType = 'okrank';
    else
        if (strfind(handles.AllFileTypes{get(hObject, 'Value')}, 'Wav'))
            handles.FileType = 'wav';
        end
    end
end
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



function OutputDirEdit_Callback(hObject, eventdata, handles)
% hObject    handle to OutputDirEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of OutputDirEdit as text
%        str2double(get(hObject,'String')) returns contents of OutputDirEdit as a double

handles.OutputDir = get(hObject, 'String');
if (~exist(handles.OutputDir, 'dir'))
    eval(['!mkdir ''' handles.OutputDir '''']);
end

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function OutputDirEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to OutputDirEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in DoAnalysisButton.
function DoAnalysisButton_Callback(hObject, eventdata, handles, varargin)
% hObject    handle to DoAnalysisButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Make directories for storing template match files
% First for normal comparisons
handles = MakeMicrolesionTemplateAnalysisDirectories(handles, 'TemplateMatchResults', 'Motif', 'Normal');
handles = MakeMicrolesionTemplateAnalysisDirectories(handles, 'TemplateMatchResults', 'Syll', 'Normal');

if (nargin <= 3)
    Temp = inputdlg('How many days do you want to analyse', 'No of days');
    handles.NoofDaysToAnalyse = str2double(Temp{1});

    FileSep = filesep;

    for i = 1:handles.NoofDaysToAnalyse,
        handles.DataDirectories{i} = uigetdir(pwd, ['Pick data directory for day #', num2str(i)]);

        [FileName, PathName] = uigetfile('*', ['Pick a filelist for all the directed songs for day #', num2str(i)]);
        handles.DirSongFileList{i} = [PathName, FileName];

        [FileName, PathName] = uigetfile('*', ['Pick a filelist for all the undirected songs for day #', num2str(i)]);
        handles.UnDirSongFileList{i} = [PathName, FileName];
    end
else
    if (strfind(varargin{1}, 'Selectivity'))
        
    else
        handles.NoofDaysToAnalyse = varargin{1};
        handles.DataDirectories = varargin{2};
        handles.DirSongFileList{i} = varargin{3};
        handles.UnDirSongFileList{i} = varargin{4};
    end
end
FileSep = filesep;
if (handles.OutputDir(end) ~= FileSep)
    handles.OutputDir(end+1) = FileSep;
end

for i = 1:handles.NoofDaysToAnalyse,
    TemplateMatch(handles.DataDirectories{i}, handles.DirSongFileList{i}, handles.FileType, handles.MotifTemplate, 'Motif', handles.DirNormal_MotifTemplateMatchOutputFilesDir, handles.ProgressText, handles.ProgressAxis, ['Day #', num2str(i), ' Dir: '], 'r', handles.TemplateTypes{handles.TemplateTypeChoiceIndex});
    TemplateMatch(handles.DataDirectories{i}, handles.UnDirSongFileList{i}, handles.FileType, handles.MotifTemplate, 'Motif', handles.UnDirNormal_MotifTemplateMatchOutputFilesDir, handles.ProgressText, handles.ProgressAxis, ['Day #', num2str(i), ' UnDir: '], 'b', handles.TemplateTypes{handles.TemplateTypeChoiceIndex});
    
    for j = 1:length(handles.SyllTemplates.SyllableTemplates),
        NoofTemplates = length(handles.SyllTemplates.SyllableTemplates{j});
        TemplateChoice = min(NoofTemplates, 2);
        TemplateMatch(handles.DataDirectories{i}, handles.DirSongFileList{i}, handles.FileType, handles.SyllTemplates.SyllableTemplates{j}{TemplateChoice}, handles.SyllTemplates.SyllableTemplates{j}{TemplateChoice}.MotifTemplate(1).Label, handles.DirNormal_SyllTemplateMatchOutputFilesDir, handles.ProgressText, handles.ProgressAxis, ['Day #', num2str(i), ' Dir: '], 'r', handles.TemplateTypes{handles.TemplateTypeChoiceIndex});
        TemplateMatch(handles.DataDirectories{i}, handles.UnDirSongFileList{i}, handles.FileType, handles.SyllTemplates.SyllableTemplates{j}{TemplateChoice}, handles.SyllTemplates.SyllableTemplates{j}{TemplateChoice}.MotifTemplate(1).Label, handles.UnDirNormal_SyllTemplateMatchOutputFilesDir, handles.ProgressText, handles.ProgressAxis, ['Day #', num2str(i), ' UnDir: '], 'b', handles.TemplateTypes{handles.TemplateTypeChoiceIndex});
    end
end

guidata(hObject, handles);
disp('Finished Analysis');


% --- Executes on button press in ResultsButton.
function ResultsButton_Callback(hObject, eventdata, handles)
% hObject    handle to ResultsButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

MicroLesionResultsAnalysis(handles);


% --- Executes on button press in SaveDataButton.
function SaveDataButton_Callback(hObject, eventdata, handles)
% hObject    handle to SaveDataButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
save([handles.OutputDir, get(handles.MotifTemplateFileNameText, 'String'), '.Results.mat'], 'handles');
disp('Finished saving data');

% --- Executes on button press in DoTemplateSelectivityAnalysis.
function DoTemplateSelectivityAnalysis_Callback(hObject, eventdata, handles)
% hObject    handle to DoTemplateSelectivityAnalysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Make directories for storing template match files
% Random song comparisons
handles = MakeMicrolesionTemplateAnalysisDirectories(handles, 'TemplateMatchResults', 'Motif', 'RandomSongComparisons');
handles = MakeMicrolesionTemplateAnalysisDirectories(handles, 'TemplateMatchResults', 'Syll', 'RandomSongComparisons');

% Make directories for storing template match files
% Random sound comparisons
handles = MakeMicrolesionTemplateAnalysisDirectories(handles, 'TemplateMatchResults', 'Motif', 'RandomSoundComparison');
handles = MakeMicrolesionTemplateAnalysisDirectories(handles, 'TemplateMatchResults', 'Syll', 'RandomSoundComparison');

% Make directories for storing template match files
% Shuffled song comparisons
handles = MakeMicrolesionTemplateAnalysisDirectories(handles, 'TemplateMatchResults', 'Motif', 'ShuffledSongComparisons');
handles = MakeMicrolesionTemplateAnalysisDirectories(handles, 'TemplateMatchResults', 'Syll', 'ShuffledSongComparisons');

% Do the random sound comparisons
FileSep = filesep;
if (handles.UnDirRandomSoundComparison_MotifTemplateMatchOutputFilesDir(end) ~= FileSep)
    handles.UnDirRandomSoundComparison_MotifTemplateMatchOutputFilesDir(end+1) = FileSep;
end

handles.RandomSoundsDir = '/home/raghav/RandomSounds/';
handles.RandomSoundFileList = '/home/raghav/RandomSounds/RandomSounds.txt';

TemplateMatch(handles.RandomSoundsDir, handles.RandomSoundFileList, 'wav', handles.MotifTemplate, 'Motif', handles.UnDirRandomSoundComparison_MotifTemplateMatchOutputFilesDir, handles.ProgressText, handles.ProgressAxis, ['Day #', num2str(i), ' UnDir: '], 'b', handles.TemplateTypes{handles.TemplateTypeChoiceIndex});
    
for j = 1:length(handles.SyllTemplates.SyllableTemplates),
    NoofTemplates = length(handles.SyllTemplates.SyllableTemplates{j});
    TemplateChoice = min(NoofTemplates, 2);
    TemplateMatch(handles.RandomSoundsDir, handles.RandomSoundFileList, 'wav', handles.SyllTemplates.SyllableTemplates{j}{TemplateChoice}, handles.SyllTemplates.SyllableTemplates{j}{TemplateChoice}.MotifTemplate(1).Label, handles.UnDirRandomSoundComparison_SyllTemplateMatchOutputFilesDir, handles.ProgressText, handles.ProgressAxis, ['Day #', num2str(i), ' Dir: '], 'r', handles.TemplateTypes{handles.TemplateTypeChoiceIndex});
end

% Do the random song comparisons
FileSep = filesep;
if (handles.UnDirRandomSongComparisons_MotifTemplateMatchOutputFilesDir(end) ~= FileSep)
    handles.UnDirRandomSongComparisons_MotifTemplateMatchOutputFilesDir(end+1) = FileSep;
end

% First for the wav files
handles.RandomSongsWavDir = '/data2/raghav/DataFromHardDrive/HVC_Microlesions/RandomComparisonSongs/WavFiles/';
handles.RandomSongsWavFileList = '/data2/raghav/DataFromHardDrive/HVC_Microlesions/RandomComparisonSongs/WavSongFiles.txt';

TemplateMatch(handles.RandomSongsWavDir, handles.RandomSongsWavFileList, 'wav', handles.MotifTemplate, 'Motif', handles.UnDirRandomSongComparisons_MotifTemplateMatchOutputFilesDir, handles.ProgressText, handles.ProgressAxis, ['Day #', num2str(i), ' UnDir: '], 'b', handles.TemplateTypes{handles.TemplateTypeChoiceIndex});
    
for j = 1:length(handles.SyllTemplates.SyllableTemplates),
    NoofTemplates = length(handles.SyllTemplates.SyllableTemplates{j});
    TemplateChoice = min(NoofTemplates, 2);
    TemplateMatch(handles.RandomSongsWavDir, handles.RandomSongsWavFileList, 'wav', handles.SyllTemplates.SyllableTemplates{j}{TemplateChoice}, handles.SyllTemplates.SyllableTemplates{j}{TemplateChoice}.MotifTemplate(1).Label, handles.UnDirRandomSongComparisons_SyllTemplateMatchOutputFilesDir, handles.ProgressText, handles.ProgressAxis, ['Day #', num2str(i), ' Dir: '], 'r', handles.TemplateTypes{handles.TemplateTypeChoiceIndex});
end

% Next for the observer files
handles.RandomSongsObsDir = '/data2/raghav/DataFromHardDrive/HVC_Microlesions/RandomComparisonSongs/ObsFiles/';
handles.RandomSongsObsFileList = '/data2/raghav/DataFromHardDrive/HVC_Microlesions/RandomComparisonSongs/ObsSongFiles.txt';

TemplateMatch(handles.RandomSongsObsDir, handles.RandomSongsObsFileList, 'obs', handles.MotifTemplate, 'Motif', handles.UnDirRandomSongComparisons_MotifTemplateMatchOutputFilesDir, handles.ProgressText, handles.ProgressAxis, ['Day #', num2str(i), ' UnDir: '], 'b', handles.TemplateTypes{handles.TemplateTypeChoiceIndex});
    
for j = 1:length(handles.SyllTemplates.SyllableTemplates),
    NoofTemplates = length(handles.SyllTemplates.SyllableTemplates{j});
    TemplateChoice = min(NoofTemplates, 2);
    TemplateMatch(handles.RandomSongsObsDir, handles.RandomSongsObsFileList, 'obs', handles.SyllTemplates.SyllableTemplates{j}{TemplateChoice}, handles.SyllTemplates.SyllableTemplates{j}{TemplateChoice}.MotifTemplate(1).Label, handles.UnDirRandomSongComparisons_SyllTemplateMatchOutputFilesDir, handles.ProgressText, handles.ProgressAxis, ['Day #', num2str(i), ' Dir: '], 'r', handles.TemplateTypes{handles.TemplateTypeChoiceIndex});
end

% Finally for the okrank files
handles.RandomSongsOKrankDir = '/data2/raghav/DataFromHardDrive/HVC_Microlesions/RandomComparisonSongs/OKrankFiles/';
handles.RandomSongsOKrankFileList = '/data2/raghav/DataFromHardDrive/HVC_Microlesions/RandomComparisonSongs/OKrankSongFiles.txt';

TemplateMatch(handles.RandomSongsOKrankDir, handles.RandomSongsOKrankFileList, 'okrank', handles.MotifTemplate, 'Motif', handles.UnDirRandomSongComparisons_MotifTemplateMatchOutputFilesDir, handles.ProgressText, handles.ProgressAxis, ['Day #', num2str(i), ' UnDir: '], 'b', handles.TemplateTypes{handles.TemplateTypeChoiceIndex});
    
for j = 1:length(handles.SyllTemplates.SyllableTemplates),
    NoofTemplates = length(handles.SyllTemplates.SyllableTemplates{j});
    TemplateChoice = min(NoofTemplates, 2);
    TemplateMatch(handles.RandomSongsOKrankDir, handles.RandomSongsOKrankFileList, 'okrank', handles.SyllTemplates.SyllableTemplates{j}{TemplateChoice}, handles.SyllTemplates.SyllableTemplates{j}{TemplateChoice}.MotifTemplate(1).Label, handles.UnDirRandomSongComparisons_SyllTemplateMatchOutputFilesDir, handles.ProgressText, handles.ProgressAxis, ['Day #', num2str(i), ' Dir: '], 'r', handles.TemplateTypes{handles.TemplateTypeChoiceIndex});
end

%=========================================

% Do the shuffled song comparisons
FileSep = filesep;
if (handles.UnDirShuffledSongComparisons_MotifTemplateMatchOutputFilesDir(end) ~= FileSep)
    handles.UnDirShuffledSongComparisons_MotifTemplateMatchOutputFilesDir(end+1) = FileSep;
end

RandomTemplateMatch(handles.DataDirectories{1}, handles.UnDirSongFileList{1}, handles.FileType, handles.MotifTemplate, 'Motif', handles.UnDirShuffledSongComparisons_MotifTemplateMatchOutputFilesDir, handles.ProgressText, handles.ProgressAxis, ['Day #', num2str(1), ' UnDir: '], 'b', handles.TemplateTypes{handles.TemplateTypeChoiceIndex});
    
for j = 1:length(handles.SyllTemplates.SyllableTemplates),
    NoofTemplates = length(handles.SyllTemplates.SyllableTemplates{j});
    TemplateChoice = min(NoofTemplates, 2);
    RandomTemplateMatch(handles.DataDirectories{1}, handles.UnDirSongFileList{1}, handles.FileType, handles.SyllTemplates.SyllableTemplates{j}{TemplateChoice}, handles.SyllTemplates.SyllableTemplates{j}{TemplateChoice}.MotifTemplate(1).Label, handles.UnDirShuffledSongComparisons_SyllTemplateMatchOutputFilesDir, handles.ProgressText, handles.ProgressAxis, ['Day #', num2str(i), ' UnDir: '], 'b', handles.TemplateTypes{handles.TemplateTypeChoiceIndex});
end

guidata(hObject, handles);
disp('Finished Analysis');

% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in TempateTypeChoiceMenu.
function TempateTypeChoiceMenu_Callback(hObject, eventdata, handles)
% hObject    handle to TempateTypeChoiceMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns TempateTypeChoiceMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from TempateTypeChoiceMenu

handles.TemplateTypes = cellstr(get(hObject, 'String'));
handles.TemplateTypeChoiceIndex = get(hObject, 'Value');
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function TempateTypeChoiceMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TempateTypeChoiceMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ExcludeBirdsEdit_Callback(hObject, eventdata, handles)
% hObject    handle to ExcludeBirdsEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ExcludeBirdsEdit as text
%        str2double(get(hObject,'String')) returns contents of ExcludeBirdsEdit as a double

Temp = get(hObject, 'String');
Commas = find(Temp == ',');
if (~isempty(Commas))
    handles.ExcludeBirds{1} = Temp(1:(Commas(1)-1));
    for i = 1:length(Commas)-1,
        handles.ExcludeBirds{i+1} = Temp(Commas(i)+1:(Commas(i+1)-1));
    end
    handles.ExcludeBirds{end+1} = Temp(Commas(end)+1:end);
else
    handles.ExcludeBirds{1} = Temp;
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function ExcludeBirdsEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ExcludeBirdsEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in LoadDetailsButton.
function LoadDetailsButton_Callback(hObject, eventdata, handles)
% hObject    handle to LoadDetailsButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[SavedFileName, SavedFileDir] = uigetfile('*.mat', 'Choose a saved data file');
load(fullfile(SavedFileDir, SavedFileName));
guidata(hObject, handles);

set(handles.BirdNameEdit, 'String', handles.BirdName);
set(handles.TempateTypeChoiceMenu, 'Value', handles.TemplateTypeChoiceIndex);

ZeroStretchIndex = find(([handles.MotifTemplate.MotifTemplate.TimeStretch] == 0) & ([handles.MotifTemplate.MotifTemplate.FreqStretch] == 0));

if (~isempty(ZeroStretchIndex))
    axes(handles.MLMATemplateAxis);
    contourf(handles.MotifTemplate.MotifTemplate(ZeroStretchIndex).MotifTemplate);
    colorbar;
    set(gca, 'YTick', []);
    set(gca, 'XTick', []);
    set(handles.MotifTemplateFileNameText, 'String', TemplateFileName);
end

if (strfind(handles.FileType, 'obs'))
    set(handles.FileTypeMenu, 'Value', 1);
else
    if (strfind(handles.FileType, 'okrank'))
        set(handles.FileTypeMenu, 'Value', 2);
    else
        if (strfind(handles.FileType, 'wav'))
            set(handles.FileTypeMenu, 'Value', 3);
        end
    end
end

if (~isempty(handles.ExcludeBirds))
    OutputString = [handles.ExcludeBirds{1}];
    for i = 2:length(handles.ExcludeBirds),
        OutputString = [OutputString, ',', handles.ExcludeBirds{i}]; 
    end
end

Index = strfind(SavedFileName, '.Results.mat');
set(handles.MotifTemplateFileNameText, 'String', SavedFileName(1:Index-1));


% --- Executes on button press in ChooseNoteFileDirButton.
function ChooseNoteFileDirButton_Callback(hObject, eventdata, handles)
% hObject    handle to ChooseNoteFileDirButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.NoteFileDir = uigetdir(pwd, 'Choose the directory with the note files');

guidata(hObject, handles);

% --- Executes on button press in ChoosePreLabelledFileListButton.
function ChoosePreLabelledFileListButton_Callback(hObject, eventdata, handles)
% hObject    handle to ChoosePreLabelledFileListButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[handles.PreLabelledFileListName, handles.PreLabelledDirName] = uigetfile('*', 'Choose pre-labelled file list');
 
FileSep = filesep;

Fid = fopen([handles.PreLabelledDirName, FileSep, handles.PreLabelledFileListName], 'r');
TempFiles = textscan(Fid, '%s', 'DeLimiter', '\n');
handles.PreLabelledFileNames = TempFiles{1};
fclose(Fid);

guidata(hObject, handles);

% --- Executes on button press in LoadPreLabelledFilesButton.
function LoadPreLabelledFilesButton_Callback(hObject, eventdata, handles)
% hObject    handle to LoadPreLabelledFilesButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
FileSep = filesep;
handles.PreLabelledTemplateMatchLabels = [];
handles.PreLabelledTemplateMatchValues = [];
handles.PreLabelledSAPFeatValues = [];

for i = 1:length(handles.PreLabelledFileNames),
    handles.Notes{i} = load([handles.NoteFileDir, FileSep, handles.PreLabelledFileNames{i}]);
    if (strfind(handles.FileType, 'okrank'))
        [Song, Fs] = ReadOKrankData(handles.PreLabelledRawDataDir, handles.PreLabelledFileNames{i}(1:end-8), 1);
        Song = Song/10;
    else
        if (strfind(handles.FileType, 'wav'))
            PresentDir = pwd;
            cd(handles.PreLabelledRawDataDir);
            [Song, Fs] = wavread(SongFile);
            cd(PresentDir);
        else
            if (strfind(handles.FileType, 'obs'))
                channel_string = strcat('obs',num2str(0),'r');
                [Song, Fs] = soundin_copy([handles.PreLabelledRawDataDir, FileSep], handles.PreLabelledFileNames{i}(1:end-8), channel_string);
                % Convert to uV - 5V on the data acquisition is 32768
                Song = Song * 5/32768;
                Song = Song/5;
            end
        end
    end
    [Feats] = CalculateSAPFeatsWithOnsets(Song, (1:1:length(Song))/Fs, Fs, handles.Notes{i}.onsets/1000, handles.Notes{i}.offsets/1000);
    Feats(:,2) = [];
    
    for j = 1:length(handles.Notes{i}.labels),
        TemplateMatchFile = dir([handles.UnDirNormal_SyllTemplateMatchOutputFilesDir, FileSep, handles.PreLabelledFileNames{i}(1:end-8), '.', handles.Notes{i}.labels(j), '.TempMatch.mat']);
        if (~isempty(TemplateMatchFile))
            PresentDir = pwd;
            cd(handles.UnDirNormal_SyllTemplateMatchOutputFilesDir);
            TemplateMatch = load(TemplateMatchFile(1).name);
            cd(PresentDir);
            TimeIndices = find((TemplateMatch.Bout.T >= (handles.Notes{i}.onsets(j)/1000 - 0.01)) & (TemplateMatch.Bout.T <= (handles.Notes{i}.onsets(j)/1000 + 0.01)));
            if (TimeIndices(1) < length(TemplateMatch.Bout.MaxBoutSeqMatch))
                handles.PreLabelledTemplateMatchLabels(end+1) = handles.Notes{i}.labels(j);
                handles.PreLabelledSAPFeatValues(end+1,:) = Feats(j,:);
                
                if (TimeIndices(end) > length(TemplateMatch.Bout.MaxBoutSeqMatch))
                    handles.PreLabelledTemplateMatchValues(end+1) = max(TemplateMatch.Bout.MaxBoutSeqMatch(TimeIndices(1):end));
                else
                    handles.PreLabelledTemplateMatchValues(end+1) = max(TemplateMatch.Bout.MaxBoutSeqMatch(TimeIndices));
                end
            end
        else
            handles.PreLabelledTemplateMatchLabels(end+1) = handles.Notes{i}.labels(j);
            handles.PreLabelledSAPFeatValues(end+1,:) = Feats(j,:);
        end
    end
end

handles.PreLabelledUniqueSylls = unique(char(handles.PreLabelledTemplateMatchLabels(find((handles.PreLabelledTemplateMatchLabels >= 97) & (handles.PreLabelledTemplateMatchLabels <= 122)))));

for i = 1:length(handles.PreLabelledUniqueSylls),
    MatchSylls = find(char(handles.PreLabelledTemplateMatchLabels) == handles.PreLabelledUniqueSylls(i));
    if (length(MatchSylls) >= 50)
        ValidSyll(i) = 1;
    else
        ValidSyll(i) = 0;
    end
end
handles.PreLabelledUniqueSylls = handles.PreLabelledUniqueSylls(find(ValidSyll == 1));

if (isfield(handles, 'MahalDist'))
    handles = rmfield(handles, 'MahalDist');
end

for i = 1:length(handles.PreLabelledUniqueSylls),
    MatchSylls = find(char(handles.PreLabelledTemplateMatchLabels) == handles.PreLabelledUniqueSylls(i));
    if (~isempty(handles.PreLabelledTemplateMatchValues))
        handles.TemplateMatchSyllThresholds{i} = mean(handles.PreLabelledTemplateMatchValues(MatchSylls)) - 2*std(handles.PreLabelledTemplateMatchValues(MatchSylls));
    end
    for j = 1:length(MatchSylls),
        handles.MahalDist{i}(j) = mahal(handles.PreLabelledSAPFeatValues(MatchSylls(j),:), handles.PreLabelledSAPFeatValues(setdiff(MatchSylls, j, 'stable'), :));
    end
end

guidata(hObject, handles);


% --- Executes on button press in ChoosePreLabelledFilesRawDataDirButton.
function ChoosePreLabelledFilesRawDataDirButton_Callback(hObject, eventdata, handles)
% hObject    handle to ChoosePreLabelledFilesRawDataDirButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.PreLabelledRawDataDir = uigetdir(pwd, 'Choose the raw data directory for the pre-labelled files');

guidata(hObject, handles);


% --- Executes on button press in SegmentAronovFeeButton.
function SegmentAronovFeeButton_Callback(hObject, eventdata, handles)
% hObject    handle to SegmentAronovFeeButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

FileSep = filesep;

PresentDir = pwd;
cd(handles.OutputDir);
handles.NoteFileOutputDir = [handles.OutputDir, handles.BirdName, '.NoteFiles'];
mkdir(handles.NoteFileOutputDir);
mkdir([handles.NoteFileOutputDir, FileSep, 'Dir']);
mkdir([handles.NoteFileOutputDir, FileSep, 'UnDir']);
cd(PresentDir);

for i = 1:length(handles.UnDirSongFileList),
    Fid = fopen(handles.UnDirSongFileList{i}, 'r');
    Temp = textscan(Fid, '%s', 'delimiter', '\n');
    fclose(Fid);

    SongFileNames = Temp{1};

    for j = 1:length(SongFileNames),
        Slash = find((SongFileNames{j} == '/') | (SongFileNames{j} == '\'));
        if (~isempty(Slash))
            SongFileNames{j} = SongFileNames{j}(Slash(end)+1:end);
        end
        if (exist([handles.NoteFileOutputDir, FileSep, 'UnDir', FileSep, SongFileNames{j}, '.not.mat'], 'file'))
            continue;
        end
        
        if (strfind(handles.FileType, 'okrank'))
            [Song, Fs] = ReadOKrankData(handles.DataDirectories{i}, SongFileNames{j}, 1);
            Song = Song/10;
        else
            if (strfind(handles.FileType, 'wav'))
                PresentDir = pwd;
                cd(handles.DataDirectories{i});
                [Song, Fs] = wavread(SongFileNames{j});
                cd(PresentDir);
            else
                if (strfind(handles.FileType, 'obs'))
                    channel_string = strcat('obs',num2str(0),'r');
                    [Song, Fs] = soundin_copy([handles.DataDirectories{i}, FileSep], SongFileNames{j}, channel_string);
                    % Convert to uV - 5V on the data acquisition is 32768
                    Song = Song * 5/32768;
                    Song = Song/5;
                end
            end
        end
        Time = (1:1:length(Song))/Fs;
    
        [LogAmplitude] = ASSLCalculateLogAmplitudeAronovFee(Song, Fs, Time, 8, 0.5);
    
        handles.UnDirThreshold{i}{j} = ASSLCalculateFisherThreshold(LogAmplitude);
        [handles.UnDirOnsets{i}{j}, handles.UnDirOffsets{i}{j}] = ASSLSegmentDataAronovFee(LogAmplitude, Fs, 7, 7, handles.UnDirThreshold{i}{j});
        handles.UnDirLabels{i}{j} = repmat('0', 1, length(handles.UnDirOnsets{i}{j})); 
        onsets = handles.UnDirOnsets{i}{j};
        offsets = handles.UnDirOffsets{i}{j};
        labels = handles.UnDirLabels{i}{j};
        min_dur = 7;
        min_int = 7;
        sm_win = 2.5;
        threshold = handles.UnDirThreshold{i}{j};
        save([handles.NoteFileOutputDir, FileSep, 'UnDir', FileSep, SongFileNames{j}, '.not.mat'], 'onsets', 'offsets', 'labels', 'min_int', 'min_dur', 'sm_win', 'threshold');
    end
end

for i = 1:length(handles.DirSongFileList),
    Fid = fopen(handles.DirSongFileList{i}, 'r');
    Temp = textscan(Fid, '%s', 'delimiter', '\n');
    fclose(Fid);

    SongFileNames = Temp{1};

    for j = 1:length(SongFileNames),
        Slash = find((SongFileNames{j} == '/') | (SongFileNames{j} == '\'));
        if (~isempty(Slash))
            SongFileNames{j} = SongFileNames{j}(Slash(end)+1:end);
        end
        if (exist([handles.NoteFileOutputDir, FileSep, 'Dir', FileSep, SongFileNames{j}, '.not.mat'], 'file'))
            continue;
        end
        if (strfind(handles.FileType, 'okrank'))
            [Song, Fs] = ReadOKrankData(handles.DataDirectories{i}, SongFileNames{j}, 1);
            Song = Song/10;
        else
            if (strfind(handles.FileType, 'wav'))
                PresentDir = pwd;
                cd(handles.DataDirectories{i});
                [Song, Fs] = wavread(SongFileNames{j});
                cd(PresentDir);
            else
                if (strfind(handles.FileType, 'obs'))
                    channel_string = strcat('obs',num2str(0),'r');
                    [Song, Fs] = soundin_copy([handles.DataDirectories{i}, FileSep], SongFileNames{j}, channel_string);
                    % Convert to uV - 5V on the data acquisition is 32768
                    Song = Song * 5/32768;
                    Song = Song/5;
                end
            end
        end
        Time = (1:1:length(Song))/Fs;
    
        [LogAmplitude] = ASSLCalculateLogAmplitudeAronovFee(Song, Fs, Time, 8, 0.5);
    
        handles.DirThreshold{i}{j} = ASSLCalculateFisherThreshold(LogAmplitude);
        [handles.DirOnsets{i}{j}, handles.DirOffsets{i}{j}] = ASSLSegmentDataAronovFee(LogAmplitude, Fs, 7, 7, handles.DirThreshold{i}{j});
        handles.DirLabels{i}{j} = repmat('0', 1, length(handles.DirOnsets{i}{j})); 
        onsets = handles.DirOnsets{i}{j};
        offsets = handles.DirOffsets{i}{j};
        labels = handles.DirLabels{i}{j};
        min_dur = 7;
        min_int = 7;
        sm_win = 2.5;
        threshold = handles.DirThreshold{i}{j};
        save([handles.NoteFileOutputDir, FileSep, 'Dir', FileSep, SongFileNames{j}, '.not.mat'], 'onsets', 'offsets', 'labels', 'min_int', 'min_dur', 'sm_win', 'threshold');
    end
end


% --- Executes on button press in ComparisonPreLabelledResultsButton.
function ComparisonPreLabelledResultsButton_Callback(hObject, eventdata, handles)
% hObject    handle to ComparisonPreLabelledResultsButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
for i = 1:length(handles.UnDirSongFileList),
    Fid = fopen(handles.UnDirSongFileList{i}, 'r');
    Temp = textscan(Fid, '%s', 'delimiter', '\n');
    fclose(Fid);

    SongFileNames = Temp{1};

    for j = 1:length(SongFileNames),
        Slash = find((SongFileNames{j} == '/') | (SongFileNames{j} == '\'));
        if (~isempty(Slash))
            SongFileNames{j} = SongFileNames{j}(Slash(end)+1:end);
        end
        if (exist([handles.NoteFileOutputDir, FileSep, 'UnDir', FileSep, SongFileNames{j}, '.not.mat'], 'file'))
            continue;
        end
        
        if (strfind(handles.FileType, 'okrank'))
            [Song, Fs] = ReadOKrankData(handles.DataDirectories{i}, SongFileNames{j}, 1);
            Song = Song/10;
        else
            if (strfind(handles.FileType, 'wav'))
                PresentDir = pwd;
                cd(handles.DataDirectories{i});
                [Song, Fs] = wavread(SongFileNames{j});
                cd(PresentDir);
            else
                if (strfind(handles.FileType, 'obs'))
                    channel_string = strcat('obs',num2str(0),'r');
                    [Song, Fs] = soundin_copy([handles.DataDirectories{i}, FileSep], SongFileNames{j}, channel_string);
                    % Convert to uV - 5V on the data acquisition is 32768
                    Song = Song * 5/32768;
                    Song = Song/5;
                end
            end
        end
        Time = (1:1:length(Song))/Fs;
    
        [LogAmplitude] = ASSLCalculateLogAmplitudeAronovFee(Song, Fs, Time, 8, 0.5);
    
        handles.UnDirThreshold{i}{j} = ASSLCalculateFisherThreshold(LogAmplitude);
        [handles.UnDirOnsets{i}{j}, handles.UnDirOffsets{i}{j}] = ASSLSegmentDataAronovFee(LogAmplitude, Fs, 7, 7, handles.UnDirThreshold{i}{j});
        [Feats] = CalculateSAPFeatsWithOnsets(Song, (1:1:length(Song))/Fs, Fs, handles.UnDirOnsets{i}{j}/1000, handles.UnDirOffsets{i}{j}/1000);
        Feats(:,2) = [];
        handles.UnDirLabels{i}{j} = repmat('0', 1, length(handles.UnDirOnsets{i}{j})); 
        for k = 1:length(handles.PreLabelledUniqueSylls),
            MatchSylls = find(char(handles.PreLabelledTemplateMatchLabels) == handles.PreLabelledUniqueSylls(k));
            handles.UnDirMahalDist{i}{j}{k} = mahal(Feats, handles.PreLabelledSAPFeatValues(MatchSylls,:));
            Matches = find(handles.UnDirMahalDist{i}{j}{k} <= (mean(handles.MahalDist{k}) + 3*std(handles.MahalDist{k})));
            if (~isempty(Matches))
                handles.UnDirLabels{i}{j}(Matches) = handles.PreLabelledUniqueSylls(k);
            end
        end
        onsets = handles.UnDirOnsets{i}{j};
        offsets = handles.UnDirOffsets{i}{j};
        labels = handles.UnDirLabels{i}{j};
        min_dur = 7;
        min_int = 7;
        sm_win = 2.5;
        threshold = handles.UnDirThreshold{i}{j};
        save([handles.NoteFileOutputDir, FileSep, 'UnDir', FileSep, SongFileNames{j}, '.not.mat'], 'onsets', 'offsets', 'labels', 'min_int', 'min_dur', 'sm_win', 'threshold');
    end
end

for i = 1:length(handles.DirSongFileList),
    Fid = fopen(handles.DirSongFileList{i}, 'r');
    Temp = textscan(Fid, '%s', 'delimiter', '\n');
    fclose(Fid);

    SongFileNames = Temp{1};

    for j = 1:length(SongFileNames),
        Slash = find((SongFileNames{j} == '/') | (SongFileNames{j} == '\'));
        if (~isempty(Slash))
            SongFileNames{j} = SongFileNames{j}(Slash(end)+1:end);
        end
        if (exist([handles.NoteFileOutputDir, FileSep, 'Dir', FileSep, SongFileNames{j}, '.not.mat'], 'file'))
            continue;
        end
        if (strfind(handles.FileType, 'okrank'))
            [Song, Fs] = ReadOKrankData(handles.DataDirectories{i}, SongFileNames{j}, 1);
            Song = Song/10;
        else
            if (strfind(handles.FileType, 'wav'))
                PresentDir = pwd;
                cd(handles.DataDirectories{i});
                [Song, Fs] = wavread(SongFileNames{j});
                cd(PresentDir);
            else
                if (strfind(handles.FileType, 'obs'))
                    channel_string = strcat('obs',num2str(0),'r');
                    [Song, Fs] = soundin_copy([handles.DataDirectories{i}, FileSep], SongFileNames{j}, channel_string);
                    % Convert to uV - 5V on the data acquisition is 32768
                    Song = Song * 5/32768;
                    Song = Song/5;
                end
            end
        end
        Time = (1:1:length(Song))/Fs;
    
        [LogAmplitude] = ASSLCalculateLogAmplitudeAronovFee(Song, Fs, Time, 8, 0.5);
    
        handles.DirThreshold{i}{j} = ASSLCalculateFisherThreshold(LogAmplitude);
        [handles.DirOnsets{i}{j}, handles.DirOffsets{i}{j}] = ASSLSegmentDataAronovFee(LogAmplitude, Fs, 7, 7, handles.DirThreshold{i}{j});
        [Feats] = CalculateSAPFeatsWithOnsets(Song, (1:1:length(Song))/Fs, Fs, handles.DirOnsets{i}{j}/1000, handles.DirOffsets{i}{j}/1000);
        Feats(:,2) = [];
        handles.DirLabels{i}{j} = repmat('0', 1, length(handles.DirOnsets{i}{j})); 
        for k = 1:length(handles.PreLabelledUniqueSylls),
            MatchSylls = find(char(handles.PreLabelledTemplateMatchLabels) == handles.PreLabelledUniqueSylls(k));
            handles.DirMahalDist{i}{j}{k} = mahal(Feats, handles.PreLabelledSAPFeatValues(MatchSylls,:));
            Matches = find(handles.DirMahalDist{i}{j}{k} <= (mean(handles.MahalDist{k}) + 3*std(handles.MahalDist{k})));
            if (~isempty(Matches))
                handles.DirLabels{i}{j}(Matches) = handles.PreLabelledUniqueSylls(k);
            end
        end
        onsets = handles.DirOnsets{i}{j};
        offsets = handles.DirOffsets{i}{j};
        labels = handles.DirLabels{i}{j};
        min_dur = 7;
        min_int = 7;
        sm_win = 2.5;
        threshold = handles.DirThreshold{i}{j};
        save([handles.NoteFileOutputDir, FileSep, 'Dir', FileSep, SongFileNames{j}, '.not.mat'], 'onsets', 'offsets', 'labels', 'min_int', 'min_dur', 'sm_win', 'threshold');
    end
end

% --- Executes on button press in pushbutton19.
function pushbutton19_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in CorrBoutSpectralContentButton.
function CorrBoutSpectralContentButton_Callback(hObject, eventdata, handles)
% hObject    handle to CorrBoutSpectralContentButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

guidata(hObject, handles);