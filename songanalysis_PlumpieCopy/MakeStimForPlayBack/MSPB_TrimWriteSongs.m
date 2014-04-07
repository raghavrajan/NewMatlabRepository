function varargout = MSPB_TrimWriteSongs(varargin)
% MSPB_TRIMWRITESONGS MATLAB code for MSPB_TrimWriteSongs.fig
%      MSPB_TRIMWRITESONGS, by itself, creates a new MSPB_TRIMWRITESONGS or raises the existing
%      singleton*.
%
%      H = MSPB_TRIMWRITESONGS returns the handle to a new MSPB_TRIMWRITESONGS or the handle to
%      the existing singleton*.
%
%      MSPB_TRIMWRITESONGS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MSPB_TRIMWRITESONGS.M with the given input arguments.
%
%      MSPB_TRIMWRITESONGS('Property','Value',...) creates a new MSPB_TRIMWRITESONGS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MSPB_TrimWriteSongs_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MSPB_TrimWriteSongs_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MSPB_TrimWriteSongs

% Last Modified by GUIDE v2.5 07-Apr-2014 11:57:26

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MSPB_TrimWriteSongs_OpeningFcn, ...
                   'gui_OutputFcn',  @MSPB_TrimWriteSongs_OutputFcn, ...
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


% --- Executes just before MSPB_TrimWriteSongs is made visible.
function MSPB_TrimWriteSongs_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MSPB_TrimWriteSongs (see VARARGIN)

% Choose default command line output for MSPB_TrimWriteSongs
handles.output = hObject;

handles.TWS.SongSilenceChoices = [{'Song'}; {'Silence'}];
handles.TWS.SongSilenceChoice = 1;
set(handles.SongSilenceChoiceMenu, 'String', handles.TWS.SongSilenceChoices);
set(handles.SongSilenceChoiceMenu, 'Value', handles.TWS.SongSilenceChoice);

handles.TWS.DirectoryName = varargin{1};
handles.TWS.FileName = varargin{2};

axes(handles.SpectrogramAxes);
PlotSpectrogramInAxis(handles.TWS.DirectoryName, handles.TWS.FileName, 'wav', handles.SpectrogramAxes);
Temp = axis;
set(gca, 'XTick', (0:1:Temp(2)));
xlabel('Time (sec)', 'FontSize', 12);
hold on;
Temp = axis;

handles.TWS.LeftBoundary = 0;
handles.TWS.RightBoundary = Temp(2);

handles.TWS.LeftBoundaryLine = plot([handles.TWS.LeftBoundary handles.TWS.LeftBoundary], [300 8000], 'r--', 'LineWidth', 2);
handles.TWS.RightBoundaryLine = plot([handles.TWS.RightBoundary handles.TWS.RightBoundary], [300 8000], 'g--', 'LineWidth', 2);

Flag = 1;

while (Flag)
    [x, y, button] = ginput(1);
    if (isempty(x) && isempty(y) && isempty(button))
        Flag = 0;
        break;
    else
        if (button == 1)
            handles.TWS.LeftBoundary = x;
            delete(handles.TWS.LeftBoundaryLine);
            handles.TWS.LeftBoundaryLine = plot([handles.TWS.LeftBoundary handles.TWS.LeftBoundary], [300 8000], 'r--', 'LineWidth', 2);
        else
            if (button == 3)
                handles.TWS.RightBoundary = x;
                delete(handles.TWS.RightBoundaryLine);
                handles.TWS.RightBoundaryLine = plot([handles.TWS.RightBoundary handles.TWS.RightBoundary], [300 8000], 'g--', 'LineWidth', 2);
            end
        end
    end
end

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MSPB_TrimWriteSongs wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = MSPB_TrimWriteSongs_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in SongSilenceChoiceMenu.
function SongSilenceChoiceMenu_Callback(hObject, eventdata, handles)
% hObject    handle to SongSilenceChoiceMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns SongSilenceChoiceMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from SongSilenceChoiceMenu
handles.TWS.SongSilenceChoice = get(hObject, 'Value');

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function SongSilenceChoiceMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SongSilenceChoiceMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in TrimWriteButton.
function TrimWriteButton_Callback(hObject, eventdata, handles)
% hObject    handle to TrimWriteButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

FileSep = filesep;

[RawData, Fs] = ASSLGetRawData(handles.TWS.DirectoryName, handles.TWS.FileName, 'wav', 1);
DataToWrite = RawData(round(Fs*handles.TWS.LeftBoundary):round(Fs*handles.TWS.RightBoundary));
OutputDir = uigetdir(pwd, 'Choose the output directory for output files');
if (OutputDir(end) ~= FileSep)
    OutputDir(end+1) = FileSep;
end

OutputFileName = [OutputDir, handles.TWS.FileName, '.Left_', num2str(handles.TWS.LeftBoundary), 's.Right_', num2str(handles.TWS.RightBoundary), 's.', handles.TWS.SongSilenceChoices{handles.TWS.SongSilenceChoice}, '.wav'];
wavwrite(DataToWrite, Fs, 16, OutputFileName);
disp(['Finished writing data to ', OutputFileName]);


% --- Executes during object creation, after setting all properties.
function TrimWriteButton_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TrimWriteButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
