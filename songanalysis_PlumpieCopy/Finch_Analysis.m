function varargout = Finch_Analysis(varargin)
% FINCH_ANALYSIS M-file for Finch_Analysis.fig
%      FINCH_ANALYSIS, by itself, creates a new FINCH_ANALYSIS or raises the existingl
%      singleton*.
%
%      H = FINCH_ANALYSIS returns the handle to a new FINCH_ANALYSIS or the handle to
%      the existing singleton*.
%
%      FINCH_ANALYSIS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FINCH_ANALYSIS.M with the given input arguments.
%
%      FINCH_ANALYSIS('Property','Value',...) creates a new FINCH_ANALYSIS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Finch_Analysis_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Finch_Analysis_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Finch_Analysis

% Last Modified by GUIDE v2.5 06-Jan-2016 04:11:03

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Finch_Analysis_OpeningFcn, ...
                   'gui_OutputFcn',  @Finch_Analysis_OutputFcn, ...
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


% --- Executes just before Finch_Analysis is made visible.
function Finch_Analysis_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Finch_Analysis (see VARARGIN)

% Choose default command line output for Finch_Analysis
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Finch_Analysis wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Finch_Analysis_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in load.
function load_Callback(hObject, eventdata, handles)
% hObject    handle to load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.FolderName = uigetdir(); 
guidata(hObject,handles)

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[handles.FileName,handles.FileDir]=uigetfile('*.mpg','Pick an mpg file to analyze')
handles.finchobj= VideoReader(handles.FileName);
handles.TotalNoFrames=round(handles.finchobj.FrameRate)*round(handles.finchobj.Duration);
if isempty(handles.TotalNoFrames)==1
    errordlg('No frames in this video. Please choose a different video')
else
end
handles.minvalue=1;   
handles.LastImage=handles.TotalNoFrames;
handles.maxvalue=handles.LastImage;
handles.Point1=zeros(handles.LastImage,2);
handles.Point2=zeros(handles.LastImage,2);
handles.Point3=zeros(handles.LastImage,2);
handles.Point4=zeros(handles.LastImage,2);
handles.Point5=zeros(handles.LastImage,2);
i=1;
k=1;
ImageDetails= read(handles.finchobj,i);
image(ImageDetails)
handles.IndexImage=k;
handles.CurrentImage=i;
LaImg=handles.LastImage;
     img=num2str(LaImg);
set(handles.edit1,'string',[img]);
guidata(hObject,handles)


% --- Executes on button press in Analyze.
function Analyze_Callback(hObject, eventdata, handles)
% hObject    handle to Analyze (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
i=handles.CurrentImage;
k=handles.CurrentImage;
ImageDetails= read(handles.finchobj,i);
image(ImageDetails)

[Point1(k,1),Point1(k,2)]=ginput(1);
[Point2(k,1),Point2(k,2)]=ginput(1);
[Point3(k,1),Point3(k,2)]=ginput(1);
[Point4(k,1),Point4(k,2)]=ginput(1);
[Point5(k,1),Point5(k,2)]=ginput(1);
handles.Point1(k,:)=Point1(k,:);
handles.Point2(k,:)=Point2(k,:);
handles.Point3(k,:)=Point3(k,:);
handles.Point4(k,:)=Point4(k,:);
handles.Point5(k,:)=Point5(k,:);

hold on
line([handles.Point1(k,1) handles.Point2(k,1)],[handles.Point1(k,2) handles.Point2(k,2)],'Marker','*','LineStyle','-')
hold on
line ([handles.Point2(k,1) handles.Point3(k,1)],[handles.Point2(k,2) handles.Point3(k,2)],'Marker','*','LineStyle','-')
hold on
line ([handles.Point3(k,1) handles.Point4(k,1)],[handles.Point3(k,2) handles.Point4(k,2)],'Marker','*','LineStyle','-')
hold on
line ([handles.Point4(k,1) handles.Point5(k,1)],[handles.Point4(k,2) handles.Point5(k,2)],'Marker','*','LineStyle','-')

guidata(hObject,handles)

% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Point1=handles.Point1;
Point2=handles.Point2;
Point3=handles.Point3;
Point4=handles.Point4;
Point5=handles.Point5;
FileName=inputdlg('Enter the name of file to be saved','Save As');
FileType='.mat';
File=horzcat(char(FileName),char(FileType));
save (File, 'Point1', 'Point2', 'Point3','Point4','Point5')

% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
i=handles.CurrentImage;
if i>=handles.LastImage
    errordlg('This is the Last image. You cannot go forward')
else
i=i+1;
img=num2str(i);
set(handles.StartImage,'string',[img]);
k=handles.CurrentImage;
k=k+1;
handles.CurrentImage=i;
handles.ImageIndex=k;
ImageDetails= read(handles.finchobj,i);
image(ImageDetails)

if handles.Point1(k,1)~=0
hold on
line([handles.Point1(k,1) handles.Point2(k,1)],[handles.Point1(i,2) handles.Point2(k,2)],'Marker','*','LineStyle','-')
hold on
line ([handles.Point2(k,1) handles.Point3(k,1)],[handles.Point2(i,2) handles.Point3(k,2)],'Marker','*','LineStyle','-')
hold on
line ([handles.Point3(k,1) handles.Point4(k,1)],[handles.Point3(i,2) handles.Point4(k,2)],'Marker','*','LineStyle','-')
hold on
line ([handles.Point4(k,1) handles.Point5(k,1)],[handles.Point4(i,2) handles.Point5(k,2)],'Marker','*','LineStyle','-')
else
end
end

guidata(hObject,handles)


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
i=handles.CurrentImage;
k=handles.ImageIndex;
if i==1
    errordlg('This is the first image. You cannot go back')
else
i=i-1;
img=num2str(i);
set(handles.StartImage,'string',[img]);
k=k-1;
handles.CurrentImage=i;
handles.ImageIndex=k;
ImageDetails=read(handles.finchobj,i);
image(ImageDetails)
hold on
line([handles.Point1(k,1) handles.Point2(k,1)],[handles.Point1(k,2) handles.Point2(k,2)],'Marker','*','LineStyle','-')
hold on
line ([handles.Point2(k,1) handles.Point3(k,1)],[handles.Point2(k,2) handles.Point3(k,2)],'Marker','*','LineStyle','-')
hold on
line ([handles.Point3(k,1) handles.Point4(k,1)],[handles.Point3(k,2) handles.Point4(k,2)],'Marker','*','LineStyle','-')
hold on
line ([handles.Point4(k,1) handles.Point5(k,1)],[handles.Point4(k,2) handles.Point5(k,2)],'Marker','*','LineStyle','-')

end

guidata(hObject,handles)




function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function StartImage_Callback(hObject, eventdata, handles)
% hObject    handle to StartImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of StartImage as text
%        str2double(get(hObject,'String')) returns contents of StartImage as a double
val=str2double(get(hObject,'String'))
if val>handles.LastImage
     errordlg('This is the Last image. You cannot go forward')
else
handles.CurrentImage=val;

i=handles.CurrentImage;
k=handles.CurrentImage;
handles.ImageIndex=k;
ImageDetails= read(handles.finchobj,i);
image(ImageDetails)
if handles.Point1(k,1)~=0
hold on
line([handles.Point1(k,1) handles.Point2(k,1)],[handles.Point1(k,2) handles.Point2(k,2)],'Marker','*','LineStyle','-')
hold on
line ([handles.Point2(k,1) handles.Point3(k,1)],[handles.Point2(k,2) handles.Point3(k,2)],'Marker','*','LineStyle','-')
hold on
line ([handles.Point3(k,1) handles.Point4(k,1)],[handles.Point3(k,2) handles.Point4(k,2)],'Marker','*','LineStyle','-')
hold on
line ([handles.Point4(k,1) handles.Point5(k,1)],[handles.Hand(k,2) handles.Point5(k,2)],'Marker','*','LineStyle','-')
else 
end
end
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function StartImage_CreateFcn(hObject, eventdata, handles)
% hObject    handle to StartImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
CurrImg=1;
     img=num2str(CurrImg);
set(hObject,'string',[img]);
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[handles.FileName,handles.FileDir]=uigetfile('*.mpg','Pick an mpg file to analyze')
handles.finchobj= VideoReader(handles.FileName);
handles.TotalNoFrames=handles.finchobj.NumberOfFrames;
if isempty(handles.TotalNoFrames)==1
    errordlg('No frames in this video. Please choose a different video')
else
end
handles.dataName=uigetfile('*.mat','Pick a matlab data file with the markings');
load(handles.dataName);
handles.Point1=Point1;
handles.Point2=Point2;
handles.Point3=Point3;
handles.Point4=Point4; 
handles.Point5=Point5;
handles.LastImage=handles.TotalNoFrames;
i=1;
k=1;
handles.IndexImage=k;
handles.CurrentImage=i;
ImageDetails=read(handles.finchobj,i);
image(ImageDetails)
if handles.Point1(k,1)~=0
hold on
line([handles.Point1(k,1) handles.Point2(k,1)],[handles.Point1(k,2) handles.Point2(k,2)],'Marker','*','LineStyle','-')
hold on
line ([handles.Point2(k,1) handles.Point3(k,1)],[handles.Point2(k,2) handles.Point3(k,2)],'Marker','*','LineStyle','-')
hold on
line ([handles.Point3(k,1) handles.Point4(k,1)],[handles.Point3(k,2) handles.Point4(k,2)],'Marker','*','LineStyle','-')
hold on
line ([handles.Point4(k,1) handles.Point5(k,1)],[handles.Point4(k,2) handles.Point5(k,2)],'Marker','*','LineStyle','-')
else 
end
LaImg=handles.LastImage;
     img=num2str(LaImg);
set(handles.edit1,'string',[img]);
img=num2str(i);
set(handles.StartImage,'string',[img]);
guidata(hObject,handles)
