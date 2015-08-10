function varargout = ASSLChooseSyllFFBoundaries(varargin)
% ASSLCHOOSESYLLFFBOUNDARIES MATLAB code for ASSLChooseSyllFFBoundaries.fig
%      ASSLCHOOSESYLLFFBOUNDARIES, by itself, creates a new ASSLCHOOSESYLLFFBOUNDARIES or raises the existing
%      singleton*.
%
%      H = ASSLCHOOSESYLLFFBOUNDARIES returns the handle to a new ASSLCHOOSESYLLFFBOUNDARIES or the handle to
%      the existing singleton*.
%
%      ASSLCHOOSESYLLFFBOUNDARIES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ASSLCHOOSESYLLFFBOUNDARIES.M with the given input arguments.
%
%      ASSLCHOOSESYLLFFBOUNDARIES('Property','Value',...) creates a new ASSLCHOOSESYLLFFBOUNDARIES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ASSLChooseSyllFFBoundaries_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ASSLChooseSyllFFBoundaries_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ASSLChooseSyllFFBoundaries

% Last Modified by GUIDE v2.5 10-Aug-2015 22:23:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ASSLChooseSyllFFBoundaries_OpeningFcn, ...
                   'gui_OutputFcn',  @ASSLChooseSyllFFBoundaries_OutputFcn, ...
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


% --- Executes just before ASSLChooseSyllFFBoundaries is made visible.
function ASSLChooseSyllFFBoundaries_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ASSLChooseSyllFFBoundaries (see VARARGIN)

% Choose default command line output for ASSLChooseSyllFFBoundaries
handles.output = hObject;

if (nargin > 0)
    handles.DataStruct = varargin{1};
end

handles.ASSLCSFFB.BoundaryChoices = [{'Specify limits by percent from beginning'}; {'Specify limits by percent from end'}; {'Specify limits by time from beginning (ms)'}; {'Specify limits by time from end (ms)'}];
set(handles.BoundaryChoiceMenu, 'String', handles.ASSLCSFFB.BoundaryChoices);
handles.ASSLCSFFB.BoundaryChoice = 1;
set(handles.BoundaryChoiceMenu, 'Value', handles.ASSLCSFFB.BoundaryChoice);

handles.ASSLCSFFB.StartLimit = 0;
set(handles.StartLimitEdit, 'String', num2str(handles.ASSLCSFFB.StartLimit));

handles.ASSLCSFFB.EndLimit = 100;
set(handles.EndLimitEdit, 'String', num2str(handles.ASSLCSFFB.EndLimit));

handles.ASSLCSFFB.UniqueSyllLabels = unique(handles.DataStruct.SyllIndexLabels);
handles.ASSLCSFFB.SyllIndex = 1;

handles.ASSLCSFFB.NumExamples = 3;
set(handles.NumExamplesEdit, 'String', num2str(handles.ASSLCSFFB.NumExamples));

for i = 1:length(handles.ASSLCSFFB.UniqueSyllLabels),
    handles.ASSLCSFFB.MaxFF{i} = [];
end
set(handles.MaxFFEdit, 'String', num2str(handles.ASSLCSFFB.MaxFF{handles.ASSLCSFFB.SyllIndex}));

handles.ASSLCSFFB.CurrentStartIndex = 1;

set(handles.SyllableIdentityLabel, 'String', ['Syllable ', handles.ASSLCSFFB.UniqueSyllLabels(handles.ASSLCSFFB.SyllIndex), ': #', num2str(handles.ASSLCSFFB.SyllIndex), ' of ', num2str(length(handles.ASSLCSFFB.UniqueSyllLabels)), ' syllables']);

[handles.ASSLCSFFB.FFSyllBoundaries, handles.DataStruct.FeatValues(:,end)]  = CalculatePlotSyllFFBoundaries(handles);
handles.ASSLCSFFB.FFBoundaryLimits(handles.ASSLCSFFB.SyllIndex, :) = [handles.ASSLCSFFB.StartLimit handles.ASSLCSFFB.EndLimit];

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ASSLChooseSyllFFBoundaries wait for user response (see UIRESUME)
% uiwait(handles.ASSLChooseFFBoundaries);


% --- Outputs from this function are returned to the command line.
function varargout = ASSLChooseSyllFFBoundaries_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in DoneButton.
function DoneButton_Callback(hObject, eventdata, handles)
% hObject    handle to DoneButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (handles.ASSLCSFFB.SyllIndex == length(handles.ASSLCSFFB.UniqueSyllLabels))
    ASSLMainWindow = findobj('Name', 'AutoSongSegmentLabel');
    if (~isempty(ASSLMainWindow))
        ASSLData = guidata(ASSLMainWindow);
        ASSLData.ASSL.FFSyllBoundaries = handles.ASSLCSFFB.FFSyllBoundaries;
        ASSLData.ASSL.FFBoundaryChoice = handles.ASSLCSFFB.BoundaryChoices{handles.ASSLCSFFB.BoundaryChoice};
        ASSLData.ASSL.FFUniqueSyllLabels = handles.ASSLCSFFB.UniqueSyllLabels;
        ASSLData.ASSL.FFBoundaryLimits = handles.ASSLCSFFB.FFBoundaryLimits;
        ASSLData.ASSL.MaxFF = handles.ASSLCSFFB.MaxFF;
        ASSLData.ASSL.FeatValues(:,end) = handles.DataStruct.FeatValues(:,end);
        ASSLData.ASSL.Raw.FundamentalFrequency = handles.DataStruct.Raw.FundamentalFrequency;
        guidata(ASSLMainWindow, ASSLData);
        msgbox('Returned syllable boundary information to Auto Song Segment Label');
    end
else
    handles.ASSLCSFFB.FFBoundaryLimits(handles.ASSLCSFFB.SyllIndex, :) = [handles.ASSLCSFFB.StartLimit handles.ASSLCSFFB.EndLimit];
    handles.ASSLCSFFB.SyllIndex = handles.ASSLCSFFB.SyllIndex + 1;
    
    handles.ASSLCSFFB.CurrentStartIndex = 1;
    
    handles.ASSLCSFFB.StartLimit = 0;
    set(handles.StartLimitEdit, 'String', num2str(handles.ASSLCSFFB.StartLimit));

    handles.ASSLCSFFB.EndLimit = 100;
    set(handles.EndLimitEdit, 'String', num2str(handles.ASSLCSFFB.EndLimit));
    
    set(handles.MaxFFEdit, 'String', num2str(handles.ASSLCSFFB.MaxFF{handles.ASSLCSFFB.SyllIndex}));
    
    set(handles.SyllableIdentityLabel, 'String', ['Syllable ', handles.ASSLCSFFB.UniqueSyllLabels(handles.ASSLCSFFB.SyllIndex), ': #', num2str(handles.ASSLCSFFB.SyllIndex), ' of ', num2str(length(handles.ASSLCSFFB.UniqueSyllLabels)), ' syllables']);
    [handles.ASSLCSFFB.FFSyllBoundaries, handles.DataStruct.FeatValues(:,end)]  = CalculatePlotSyllFFBoundaries(handles);
end
guidata(hObject, handles);

% --- Executes on selection change in BoundaryChoiceMenu.
function BoundaryChoiceMenu_Callback(hObject, eventdata, handles)
% hObject    handle to BoundaryChoiceMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns BoundaryChoiceMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from BoundaryChoiceMenu

handles.ASSLCSFFB.BoundaryChoice = get(hObject, 'Value');
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function BoundaryChoiceMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BoundaryChoiceMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function StartLimitEdit_Callback(hObject, eventdata, handles)
% hObject    handle to StartLimitEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of StartLimitEdit as text
%        str2double(get(hObject,'String')) returns contents of StartLimitEdit as a double

handles.ASSLCSFFB.StartLimit = str2double(get(hObject, 'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function StartLimitEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to StartLimitEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function EndLimitEdit_Callback(hObject, eventdata, handles)
% hObject    handle to EndLimitEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EndLimitEdit as text
%        str2double(get(hObject,'String')) returns contents of EndLimitEdit as a double

handles.ASSLCSFFB.EndLimit = str2double(get(hObject, 'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function EndLimitEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EndLimitEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ReCalculateButton.
function ReCalculateButton_Callback(hObject, eventdata, handles)
% hObject    handle to ReCalculateButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[handles.ASSLCSFFB.FFSyllBoundaries, handles.DataStruct.FeatValues(:,end), handles.DataStruct.Raw.FundamentalFrequency]  = CalculatePlotSyllFFBoundaries(handles, handles.ASSLCSFFB.MaxFF{handles.ASSLCSFFB.SyllIndex});
handles.ASSLCSFFB.FFBoundaryLimits(handles.ASSLCSFFB.SyllIndex, :) = [handles.ASSLCSFFB.StartLimit handles.ASSLCSFFB.EndLimit];
guidata(hObject, handles);

function NumExamplesEdit_Callback(hObject, eventdata, handles)
% hObject    handle to NumExamplesEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NumExamplesEdit as text
%        str2double(get(hObject,'String')) returns contents of NumExamplesEdit as a double

handles.ASSLCSFFB.NumExamples = str2double(get(hObject, 'String'));
set(handles.NextButton, 'String', ['Next ', num2str(handles.ASSLCSFFB.NumExamples)]);
set(handles.PrevButton, 'String', ['Prev ', num2str(handles.ASSLCSFFB.NumExamples)]);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function NumExamplesEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NumExamplesEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in NextButton.
function NextButton_Callback(hObject, eventdata, handles)
% hObject    handle to NextButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.ASSLCSFFB.CurrentStartIndex = handles.ASSLCSFFB.CurrentStartIndex + handles.ASSLCSFFB.NumExamples;
if (handles.ASSLCSFFB.CurrentStartIndex < 0)
    handles.ASSLCSFFB.CurrentStartIndex = 1;
end
ASSLPlotSyllFFBoundaries(handles);
guidata(hObject, handles);

% --- Executes on button press in PrevButton.
function PrevButton_Callback(hObject, eventdata, handles)
% hObject    handle to PrevButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.ASSLCSFFB.CurrentStartIndex = handles.ASSLCSFFB.CurrentStartIndex - handles.ASSLCSFFB.NumExamples;
if (handles.ASSLCSFFB.CurrentStartIndex < 0)
    handles.ASSLCSFFB.CurrentStartIndex = 1;
end
ASSLPlotSyllFFBoundaries(handles);
guidata(hObject, handles);



function MaxFFEdit_Callback(hObject, eventdata, handles)
% hObject    handle to MaxFFEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MaxFFEdit as text
%        str2double(get(hObject,'String')) returns contents of MaxFFEdit as a double

handles.ASSLCSFFB.MaxFF{handles.ASSLCSFFB.SyllIndex} = str2double(get(hObject, 'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function MaxFFEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MaxFFEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
