function varargout = ASSLPlotFeatures(varargin)
% ASSLPLOTFEATURES M-file for ASSLPlotFeatures.fig
%      ASSLPLOTFEATURES, by itself, creates a new ASSLPLOTFEATURES or raises the existing
%      singleton*.
%
%      H = ASSLPLOTFEATURES returns the handle to a new ASSLPLOTFEATURES or the handle to
%      the existing singleton*.
%
%      ASSLPLOTFEATURES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ASSLPLOTFEATURES.M with the given input arguments.
%
%      ASSLPLOTFEATURES('Property','Value',...) creates a new ASSLPLOTFEATURES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ASSLPlotFeatures_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ASSLPlotFeatures_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ASSLPlotFeatures

% Last Modified by GUIDE v2.5 09-Jul-2012 00:06:35

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ASSLPlotFeatures_OpeningFcn, ...
                   'gui_OutputFcn',  @ASSLPlotFeatures_OutputFcn, ...
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


% --- Executes just before ASSLPlotFeatures is made visible.
function ASSLPlotFeatures_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ASSLPlotFeatures (see VARARGIN)

% Choose default command line output for ASSLPlotFeatures
handles.output = hObject;

if (length(varargin) > 0)
    handles.ASSLPF.FeatVals = varargin{1};
    handles.ASSLPF.FeatNames = varargin{2};
    handles.DataStruct = varargin{3};
    handles.ASSLPF.SyllIndex = 1;
end
    
set(handles.XListBox, 'String', handles.ASSLPF.FeatNames);
set(handles.XListBox, 'Value', 1);
set(handles.YListBox, 'String', handles.ASSLPF.FeatNames);
set(handles.YListBox, 'Value', 1);

handles.ASSLPF.FeatNames = cellstr(get(handles.XListBox, 'String'));
handles.ASSLPF.XVal = get(handles.XListBox, 'Value');
handles.ASSLPF.YVal = get(handles.YListBox, 'Value');

ASSLPFPlotData(handles);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ASSLPlotFeatures wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ASSLPlotFeatures_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in XListBox.
function XListBox_Callback(hObject, eventdata, handles)
% hObject    handle to XListBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns XListBox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from XListBox
handles.ASSLPF.XVal = get(hObject, 'Value');
axes(handles.ASSLFeaturePlotAxis);
hold off;
plot(handles.ASSLPF.FeatVals(:,handles.ASSLPF.XVal), handles.ASSLPF.FeatVals(:,handles.ASSLPF.YVal), 'k+', 'MarkerSize', 2);
hold on;
plot(handles.ASSLPF.FeatVals(handles.ASSLPF.SyllIndex, handles.ASSLPF.XVal), handles.ASSLPF.FeatVals(handles.ASSLPF.SyllIndex, handles.ASSLPF.YVal), 'ro', 'MarkerSize', 6);
xlabel(handles.ASSLPF.FeatNames{handles.ASSLPF.XVal}, 'FontSize', 12, 'FontWeight', 'bold');
ylabel(handles.ASSLPF.FeatNames{handles.ASSLPF.YVal}, 'FontSize', 12, 'FontWeight', 'bold');
set(gca, 'FontSize', 12, 'FontWeight', 'bold');
axis tight;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function XListBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to XListBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in YListBox.
function YListBox_Callback(hObject, eventdata, handles)
% hObject    handle to YListBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns YListBox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from YListBox

handles.ASSLPF.YVal = get(hObject, 'Value');
axes(handles.ASSLFeaturePlotAxis);
hold off;
plot(handles.ASSLPF.FeatVals(:,handles.ASSLPF.XVal), handles.ASSLPF.FeatVals(:,handles.ASSLPF.YVal), 'k+', 'MarkerSize', 2);
hold on;
plot(handles.ASSLPF.FeatVals(handles.ASSLPF.SyllIndex, handles.ASSLPF.XVal), handles.ASSLPF.FeatVals(handles.ASSLPF.SyllIndex, handles.ASSLPF.YVal), 'ro', 'MarkerSize', 6);
axis tight;
xlabel(handles.ASSLPF.FeatNames{handles.ASSLPF.XVal}, 'FontSize', 12, 'FontWeight', 'bold');
ylabel(handles.ASSLPF.FeatNames{handles.ASSLPF.YVal}, 'FontSize', 12, 'FontWeight', 'bold');
set(gca, 'FontSize', 12, 'FontWeight', 'bold');
guidata(hObject, handles);
% --- Executes during object creation, after setting all properties.
function YListBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to YListBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in PrevSyllButton.
function PrevSyllButton_Callback(hObject, eventdata, handles)
% hObject    handle to PrevSyllButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.ASSLPF.SyllIndex = handles.ASSLPF.SyllIndex - 1;
if (handles.ASSLPF.SyllIndex < 1)
    handles.ASSLPF.SyllIndex = 1;
end
if (handles.ASSLPF.SyllIndex > size(handles.DataStruct.SyllIndices,1))
    handles.ASSLPF.SyllIndex = size(handles.DataStruct.SyllIndices,1);
end
guidata(hObject, handles);

ASSLPFPlotData(handles);

% --- Executes on button press in NextSyllButton.
function NextSyllButton_Callback(hObject, eventdata, handles)
% hObject    handle to NextSyllButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.ASSLPF.SyllIndex = handles.ASSLPF.SyllIndex + 1;
if (handles.ASSLPF.SyllIndex < 1)
    handles.ASSLPF.SyllIndex = 1;
end
if (handles.ASSLPF.SyllIndex > size(handles.DataStruct.SyllIndices,1))
    handles.ASSLPF.SyllIndex = size(handles.DataStruct.SyllIndices,1);
end
guidata(hObject, handles);

ASSLPFPlotData(handles);


% --- Executes on button press in MakeTemplateButton.
function MakeTemplateButton_Callback(hObject, eventdata, handles)
% hObject    handle to MakeTemplateButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

TimeRange(1) = handles.DataStruct.SyllOnsets{handles.DataStruct.SyllIndices(handles.ASSLPF.SyllIndex,1)}(handles.DataStruct.SyllIndices(handles.ASSLPF.SyllIndex,2))/1000 - 0.005;
TimeRange(2) = handles.DataStruct.SyllOffsets{handles.DataStruct.SyllIndices(handles.ASSLPF.SyllIndex,1)}(handles.DataStruct.SyllIndices(handles.ASSLPF.SyllIndex,2))/1000 + 0.005;

PresentDir = pwd;
cd(handles.DataStruct.DirName);
if (~exist('Templates', 'dir'))
    mkdir(handles.DataStruct.DirName, 'Templates');
end
cd(PresentDir);


if (ispc)
    handles.TemplateFileDir = [handles.DataStruct.DirName, '\Templates\'];
else
    handles.TemplateFileDir = [handles.DataStruct.DirName, '/Templates/'];
end

ASSLMakeTemplatesSpectralMatchAnalysis(handles.DataStruct.DirName, handles.DataStruct.FileName{handles.DataStruct.SyllIndices(handles.ASSLPF.SyllIndex,1)}, handles.DataStruct.FileType, handles.DataStruct.SongChanNo, TimeRange, handles.DataStruct.FFTWinSizeTempMatch, handles.DataStruct.FFTWinOverlapTempMatch, handles.TemplateFileDir, handles.TemplateFileName, handles.SyllableLabel)



function TemplateFileNameEdit_Callback(hObject, eventdata, handles)
% hObject    handle to TemplateFileNameEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TemplateFileNameEdit as text
%        str2double(get(hObject,'String')) returns contents of TemplateFileNameEdit as a double
handles.TemplateFileName = get(hObject, 'String');
guidata(hObject, handles);



function SyllableLabelEdit_Callback(hObject, eventdata, handles)
% hObject    handle to SyllableLabelEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SyllableLabelEdit as text
%        str2double(get(hObject,'String')) returns contents of SyllableLabelEdit as a double
handles.SyllableLabel = get(hObject, 'String');
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function SyllableLabelEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SyllableLabelEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in PickSyllableButton.
function PickSyllableButton_Callback(hObject, eventdata, handles)
% hObject    handle to PickSyllableButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[SyllX, SyllY] = ginput(1);
NearestSyll = knnsearch(handles.ASSLPF.FeatVals(:,[handles.ASSLPF.XVal handles.ASSLPF.YVal]), [SyllX SyllY], 'Distance', 'mahalanobis');

handles.ASSLPF.SyllIndex = handles.DataStruct.SyllIndices(NearestSyll,3);
guidata(hObject, handles);

ASSLPFPlotData(handles);
