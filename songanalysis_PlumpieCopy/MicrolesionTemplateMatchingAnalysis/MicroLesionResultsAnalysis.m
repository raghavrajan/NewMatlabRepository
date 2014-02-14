function varargout = MicroLesionResultsAnalysis(varargin)
% MICROLESIONRESULTSANALYSIS MATLAB code for MicroLesionResultsAnalysis.fig
%      MICROLESIONRESULTSANALYSIS, by itself, creates a new MICROLESIONRESULTSANALYSIS or raises the existing
%      singleton*.
%
%      H = MICROLESIONRESULTSANALYSIS returns the handle to a new MICROLESIONRESULTSANALYSIS or the handle to
%      the existing singleton*.
%
%      MICROLESIONRESULTSANALYSIS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MICROLESIONRESULTSANALYSIS.M with the given input arguments.
%
%      MICROLESIONRESULTSANALYSIS('Property','Value',...) creates a new MICROLESIONRESULTSANALYSIS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MicroLesionResultsAnalysis_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MicroLesionResultsAnalysis_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MicroLesionResultsAnalysis

% Last Modified by GUIDE v2.5 20-Jul-2013 13:35:59

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MicroLesionResultsAnalysis_OpeningFcn, ...
                   'gui_OutputFcn',  @MicroLesionResultsAnalysis_OutputFcn, ...
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


% --- Executes just before MicroLesionResultsAnalysis is made visible.
function MicroLesionResultsAnalysis_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MicroLesionResultsAnalysis (see VARARGIN)

% Choose default command line output for MicroLesionResultsAnalysis
handles.output = hObject;

if (nargin > 0)
    InputData = varargin{1};
    
end

ProgramDir = which('MicroLesionMainAnalysis');
SlashIndex = find(ProgramDir == filesep);
handles.AnalysisProgramsDir = [ProgramDir(1:SlashIndex(end)), 'TemplateMatchResultsAnalysis'];
cd(handles.AnalysisProgramsDir);

AnalysisFiles = dir('*.m');

for i = 1:length(AnalysisFiles),
    handles.AnalysisFiles{i} = AnalysisFiles(i).name(1:end-2);
    ListBoxString{i} = handles.AnalysisFiles{i};
end
set(handles.AnalysisListBox, 'String', ListBoxString);

SyllChoiceString{1} = InputData.MotifTemplate.MotifTemplate(1).Label;
handles.SyllableChoices{1} = InputData.MotifTemplate.MotifTemplate(1).Label;


for i = 1:length(InputData.SyllTemplates.SyllableTemplates),
    SyllChoiceString{i+1} = InputData.SyllTemplates.SyllableTemplates{i}{1}.MotifTemplate(1).Label;
    handles.SyllableChoices{i+1} = InputData.SyllTemplates.SyllableTemplates{i}{1}.MotifTemplate(1).Label;    
end
set(handles.SyllableChoiceListBox, 'String', SyllChoiceString);

handles.InputData = InputData;

OutputFile = [handles.InputData.OutputDir, get(handles.InputData.MotifTemplateFileNameText, 'String'), '.', get(handles.InputData.SyllTemplateFileNameText, 'String'), '.mat'];
if (exist(OutputFile, 'file'))
    Temp = load(OutputFile);
    handles.MatchPeaksData = Temp.handles.MatchPeaksData;
    clear Temp;
else
    handles.MatchPeaksData = AnalyzeTemplateMatchResultsData(InputData);
end

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MicroLesionResultsAnalysis wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = MicroLesionResultsAnalysis_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on selection change in AnalysisListBox.
function AnalysisListBox_Callback(hObject, eventdata, handles)
% hObject    handle to AnalysisListBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns AnalysisListBox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from AnalysisListBox

feval(handles.AnalysisFiles{get(hObject, 'Value')}, handles.MatchPeaksData, handles.SyllableChoices{get(handles.SyllableChoiceListBox, 'Value')}, handles.InputData, get(handles.SyllableChoiceListBox, 'String'));

% --- Executes during object creation, after setting all properties.
function AnalysisListBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AnalysisListBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in SyllableChoiceListBox.
function SyllableChoiceListBox_Callback(hObject, eventdata, handles)
% hObject    handle to SyllableChoiceListBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns SyllableChoiceListBox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from SyllableChoiceListBox


% --- Executes during object creation, after setting all properties.
function SyllableChoiceListBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SyllableChoiceListBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in SaveResultsButton.
function SaveResultsButton_Callback(hObject, eventdata, handles)
% hObject    handle to SaveResultsButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

OutputFile = [handles.InputData.OutputDir, get(handles.InputData.MotifTemplateFileNameText, 'String'), '.', get(handles.InputData.SyllTemplateFileNameText, 'String'), '.mat'];
save(OutputFile, 'handles');
disp('Finished saving results');

% --- Executes on button press in SavePlotResultsButton.
function SavePlotResultsButton_Callback(hObject, eventdata, handles)
% hObject    handle to SavePlotResultsButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
