function varargout = MakeTemplates(varargin)
% MAKETEMPLATES M-file for MakeTemplates.fig
%      MAKETEMPLATES, by itself, creates a new MAKETEMPLATES or raises the existing
%      singleton*.
%
%      H = MAKETEMPLATES returns the handle to a new MAKETEMPLATES or the handle to
%      the existing singleton*.
%
%      MAKETEMPLATES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MAKETEMPLATES.M with the given input arguments.
%
%      MAKETEMPLATES('Property','Value',...) creates a new MAKETEMPLATES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MakeTemplates_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MakeTemplates_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MakeTemplates

% Last Modified by GUIDE v2.5 29-Sep-2009 11:23:50

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MakeTemplates_OpeningFcn, ...
                   'gui_OutputFcn',  @MakeTemplates_OutputFcn, ...
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


% --- Executes just before MakeTemplates is made visible.
function MakeTemplates_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MakeTemplates (see VARARGIN)

% Choose default command line output for MakeTemplates
handles.output = hObject;

%==========================================================================
% Variables used in the program - Raghav 09.28.2009
%==========================================================================
handles.MakeTemplatesData.DataDirectory = pwd;
handles.MakeTemplatesData.FileType = get(handles.FileTypeEdit, 'String');
handles.MakeTemplatesData.TemplateFileName = '';
handles.MakeTemplatesData.FFTWindowLength = str2double(get(handles.FFTWindowLengthEdit, 'String'));
handles.MakeTemplatesData.FFTWindowOverlap = str2double(get(handles.FFTWindowOverlapEdit, 'String'));
%==========================================================================

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MakeTemplates wait for user response (see UIRESUME)
% uiwait(handles.MakeTemplatesMainWindow);


% --- Outputs from this function are returned to the command line.
function varargout = MakeTemplates_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in LoadSongFileButton.
function LoadSongFileButton_Callback(hObject, eventdata, handles)
% hObject    handle to LoadSongFileButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[handles.MakeTemplatesData.SongFileName, Pathname] = uigetfile('*', 'Pick a data file');

if (handles.MakeTemplatesData.DataDirectory(end) ~= '/')
    handles.MakeTemplatesData.DataDirectory(end + 1) = '/';
end
   
if (strfind(handles.MakeTemplatesData.FileType, 'okrank'))
        [RawData, Fs] = ReadOKrankData(handles.MakeTemplatesData.DataDirectory, handles.MakeTemplatesData.SongFileName, 0);
else
    if (strfind(handles.MakeTemplatesData.FileType, 'obs'))
        [RawData, Fs] = soundin_copy(handles.MakeTemplatesData.DataDirectory, handles.MakeTemplatesData.SongFileName, handles.MakeTemplatesData.FileType);    
        RawData = RawData/max(RawData);
        RawData = RawData - mean(RawData);
    else
        if (strfind(handles.MakeTemplatesData.FileType, 'wav'))
            [RawData, Fs] = wavread(handles.MakeTemplatesData.SongFileName);    
        end
    end
end

Time = (1:1:length(RawData))/Fs;

% Now using an 8 pole butterworth bandpass filter as default.
[b,a]=butter(8,[300*2/Fs, 10000*2/Fs]);

FiltSong=filtfilt(b, a, RawData);
  
if length(RawData) ~= length(FiltSong) 
  disp(['warning! bandpass: input and output file lengths do not match!']);
end

RawData = FiltSong;
clear FiltSong;

Freq = Fs/2*linspace(0, 1, handles.MakeTemplatesData.FFTWindowLength/2+1);
FreqRows = find((Freq >= 1700) & (Freq <= 7100));
[Y, F, T, P] = spectrogram(RawData, handles.MakeTemplatesData.FFTWindowLength, handles.MakeTemplatesData.FFTWindowOverlap, Freq, Fs, 'yaxis');

axes(handles.AmplitudePlot);
hold off;
plot(Time, RawData);
hold on;
plot(T, sum(abs(P))/max(sum(abs(P))), 'r');

axis tight;

set(handles.AmplitudePlot, 'FontSize', 12, 'FontWeight', 'bold');
%xlabel('Time (sec)', 'FontSize', 12, 'FontWeight', 'bold');
%ylabel('Amplitude', 'FontSize', 12, 'FontWeight', 'bold');

handles.MakeTemplatesData.AmpPlotAxis = axis;
handles.MakeTemplatesData.AmpPlotAxis(1:2) = [Time(1) Time(end)];

%MedianMotif.FileName = handles.MakeTemplatesData.SongFileName;
%MedianMotif.Length = Time(end) - Time(1);
%MedianMotif.StartTime = Time(1);

axes(handles.SpectrogramPlot);
hold off;
plot_motif_spectrogram(RawData, Fs, gcf, handles.SpectrogramPlot);
%PlotMotifSpectrogram(handles.MakeTemplatesData.DataDirectory, handles.MakeTemplatesData.FileType, MedianMotif, handles.MakeTemplatesMainWindow, handles.SpectrogramPlot);
hold on;
handles.MakeTemplatesData.SpecPlotAxis = handles.MakeTemplatesData.AmpPlotAxis;
handles.MakeTemplatesData.SpecPlotAxis(3:4) = [300 10000];
axis(handles.MakeTemplatesData.SpecPlotAxis);
%xlabel('Time (sec)', 'FontSize', 12, 'FontWeight', 'bold');
%ylabel('Frequency (Hz)', 'FontSize', 12, 'FontWeight', 'bold');

handles.MakeTemplatesData.LabelPlotAxis = handles.MakeTemplatesData.AmpPlotAxis;
handles.MakeTemplatesData.LabelPlotAxis(3:4) = [0 1];

axes(handles.LabelPlot);
plot(0, 0.5, 'w+');
axis(handles.MakeTemplatesData.LabelPlotAxis);
%set(handles.LabelPlot, 'Visible', 'off');

handles.MakeTemplatesData.PresentXCoordinates = [Time(1) Time(end)];
handles.MakeTemplatesData.RawData = RawData;
handles.MakeTemplatesData.Fs = Fs;
handles.MakeTemplatesData.Time = Time;
handles.MakeTemplatesData.Freq = Freq;
handles.MakeTemplatesData.FreqRows = FreqRows;

handles.MakeTemplatesData.SyllableBoundaries = [];
handles.MakeTemplatesData.SyllableLabels = [];
handles.SyllableLabels = [];
handles.AmpSyllBoundary = [];
handles.SpecSyllBoundary = [];

clear RawData Fs Y F T P b a Freq Time FreqRows;
guidata(handles.MakeTemplatesMainWindow, handles);

disp('Finished');


function FileTypeEdit_Callback(hObject, eventdata, handles)
% hObject    handle to FileTypeEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FileTypeEdit as text
%        str2double(get(hObject,'String')) returns contents of FileTypeEdit as a double

handles.MakeTemplatesData.FileType = get(hObject,'String');
guidata(handles.MakeTemplatesMainWindow, handles);
disp(['File type is ', handles.MakeTemplatesData.FileType]);

function FFTWindowLengthEdit_Callback(hObject, eventdata, handles)
% hObject    handle to FFTWindowLengthEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FFTWindowLengthEdit as text
%        str2double(get(hObject,'String')) returns contents of FFTWindowLengthEdit as a double

handles.MakeTemplatesData.FFTWindowLength = str2double(get(hObject,'String'));
guidata(handles.MakeTemplatesMainWindow, handles);
if (isfield(handles.MakeTemplatesData, 'Fs'))
    Fs = handles.MakeTemplatesData.Fs;
    set(handles.FFTWindowLengthText, 'String', [num2str(round(handles.MakeTemplatesData.FFTWindowLength * 1000/Fs)), ' ms']);
else
    set(handles.FFTWindowLengthText, 'String', [num2str(round(handles.MakeTemplatesData.FFTWindowLength * 1000/32000)), ' ms']);
end

% --- Executes during object creation, after setting all properties.
function FFTWindowLengthEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FFTWindowLengthEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function FFTWindowOverlapEdit_Callback(hObject, eventdata, handles)
% hObject    handle to FFTWindowOverlapEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FFTWindowOverlapEdit as text
%        str2double(get(hObject,'String')) returns contents of FFTWindowOverlapEdit as a double

handles.MakeTemplatesData.FFTWindowOverlap = str2double(get(hObject,'String'));
guidata(handles.MakeTemplatesMainWindow, handles);
if (isfield(handles.MakeTemplatesData, 'Fs'))
    Fs = handles.MakeTemplatesData.Fs;
    set(handles.FFTWindowOverlapText, 'String', [num2str(round(handles.MakeTemplatesData.FFTWindowOverlap * 1000/Fs)), ' ms']);
else
    set(handles.FFTWindowOverlapText, 'String', [num2str(round(handles.MakeTemplatesData.FFTWindowOverlap * 1000/32000)), ' ms']);
end


% --- Executes during object creation, after setting all properties.
function FFTWindowOverlapEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FFTWindowOverlapEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ChooseDirectoryButton.
function ChooseDirectoryButton_Callback(hObject, eventdata, handles)
% hObject    handle to ChooseDirectoryButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.MakeTemplatesData.DataDirectory = uigetdir('','Choose the data directory');
guidata(handles.MakeTemplatesMainWindow, handles);
cd(handles.MakeTemplatesData.DataDirectory);
disp(['Data directory is ', handles.MakeTemplatesData.DataDirectory]);

% --- Executes on button press in MarkSyllableButton.
function MarkSyllableButton_Callback(hObject, eventdata, handles)
% hObject    handle to MarkSyllableButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.ZoomInstructionsText, 'String', 'Left mouse button for left edge, right mouse button for right edge and then press enter (q to quit)');

axes(handles.AmplitudePlot);
TempAxis = axis;

handles.MakeTemplatesData.SyllableBoundaries(end+1,1:2) = [TempAxis(1) TempAxis(2)];
handles.MakeTemplatesData.SyllableLabels{end+1} = '0';

SyllBoundaryPatch = patch([TempAxis(1) TempAxis(1) TempAxis(2) TempAxis(2)], [TempAxis(3) TempAxis(4) TempAxis(4) TempAxis(3)], 'g');
set(SyllBoundaryPatch, 'FaceAlpha', 0.2);
LeftLine = plot([TempAxis(1) TempAxis(1)], [TempAxis(3) TempAxis(4)], '-.r');
RightLine = plot([TempAxis(2) TempAxis(2)], [TempAxis(3) TempAxis(4)], '-.g');
SyllableBoundary(1,:) = TempAxis(1:2);

Flag = 1;

while(Flag)
    [x, y, button] = ginput(1);
    
    if (((length(x) == 0) && (length(y) == 0) && (length(button) == 0)) || (button == 113))
        Flag = 0;
        if (exist('LeftLine', 'var'))
            delete(LeftLine);
        end
        if (exist('RightLine', 'var'))
            delete(RightLine);
        end
        if (exist('SyllBoundaryPatch', 'var'))
            delete(SyllBoundaryPatch);
        end
        if (button == 113)
            disp('Quit syllable boundary marking');
        else
            axes(handles.AmplitudePlot);
            SyllableBoundary(3:4) = SyllableBoundary(2);
            SyllableBoundary(2) = SyllableBoundary(1);
            handles.AmpSyllBoundary(end + 1) = plot(SyllableBoundary, [0 (0.9 * TempAxis(4)) (0.9 * TempAxis(4)) 0], 'm');

            axes(handles.SpectrogramPlot);
            TempAxis = axis;
            handles.SpecSyllBoundary(end + 1) = plot(SyllableBoundary, [0 (0.9 * TempAxis(4)) (0.9 * TempAxis(4)) 0], 'm');
            
            SyllLabel = inputdlg('Enter a letter for the syllable');
            handles.MakeTemplatesData.SyllableLabels{end} = SyllLabel{1};
            axes(handles.LabelPlot);
            handles.SyllableLabels(end + 1) = text([(SyllableBoundary(1) + SyllableBoundary(3))/2], 0.5, handles.MakeTemplatesData.SyllableLabels{end});
            set(handles.SyllableLabels(end), 'FontSize', 12, 'FontWeight', 'bold');
        end
            break;
    else
        if (button == 1)
            delete(LeftLine);
            axes(handles.AmplitudePlot);
            SyllableBoundary(1,1) = x;
            LeftLine = plot([x x], [TempAxis(3) TempAxis(4)], '-.r');
            delete(SyllBoundaryPatch);
            SyllBoundaryPatch = patch([SyllableBoundary(1,1) SyllableBoundary(1,1) SyllableBoundary(1,2) SyllableBoundary(1,2)], [TempAxis(3) TempAxis(4) TempAxis(4) TempAxis(3)], 'g');
            set(SyllBoundaryPatch, 'FaceAlpha', 0.2);
            handles.MakeTemplatesData.SyllableBoundaries(end, 1) = x;
        else
            if (button == 3)
                delete(RightLine);
                axes(handles.AmplitudePlot);
                RightLine = plot([x x], [TempAxis(3) TempAxis(4)], '-.g');
                handles.MakeTemplatesData.SyllableBoundaries(end, 2) = x;
                SyllableBoundary(1,2) = x;
                delete(SyllBoundaryPatch);
                SyllBoundaryPatch = patch([SyllableBoundary(1,1) SyllableBoundary(1,1) SyllableBoundary(1,2) SyllableBoundary(1,2)], [TempAxis(3) TempAxis(4) TempAxis(4) TempAxis(3)], 'g');                
                set(SyllBoundaryPatch, 'FaceAlpha', 0.2);
            end
        end
    end
end

set(handles.ZoomInstructionsText, 'String', 'Watch this space for instructions for zoom and marking syllable boundaries');
guidata(handles.MakeTemplatesMainWindow, handles);

% --- Executes on button press in SaveTemplatesButton.
function SaveTemplatesButton_Callback(hObject, eventdata, handles)
% hObject    handle to SaveTemplatesButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[Test, SortedIndices] = sort(handles.MakeTemplatesData.SyllableBoundaries(:,1));
handles.MakeTemplatesData.SyllableLabels = handles.MakeTemplatesData.SyllableLabels(SortedIndices);
handles.MakeTemplatesData.SyllableBoundaries = handles.MakeTemplatesData.SyllableBoundaries(SortedIndices,:);

for i = 1:length(handles.MakeTemplatesData.SyllableLabels);
    SyllableLabels(i) = handles.MakeTemplatesData.SyllableLabels{i};
end

UniqueSyllableLabels = unique(SyllableLabels);

for i = 1:length(UniqueSyllableLabels),
    Matches = find(SyllableLabels == UniqueSyllableLabels(i));
    Templates(i).Label = UniqueSyllableLabels(i);
    SyllableData = handles.MakeTemplatesData.RawData(round(handles.MakeTemplatesData.SyllableBoundaries(Matches(1),1) * handles.MakeTemplatesData.Fs):round(handles.MakeTemplatesData.SyllableBoundaries(Matches(1),2) * handles.MakeTemplatesData.Fs));
    [Y1, F1, T1, P1] = spectrogram(SyllableData, handles.MakeTemplatesData.FFTWindowLength, handles.MakeTemplatesData.FFTWindowOverlap, handles.MakeTemplatesData.FFTWindowLength, handles.MakeTemplatesData.Fs, 'yaxis');
    
    nfft = handles.MakeTemplatesData.FFTWindowLength;
    nfft = round(handles.MakeTemplatesData.Fs *2/1000);
    nfft = 2^nextpow2(nfft);
    spect_win = hanning(nfft);
    noverlap = round(0.5*length(spect_win)); %number of overlapping points       
%    now calculate spectrogram
    [spect, freq, time] = spectrogram(SyllableData, spect_win, noverlap, nfft, handles.MakeTemplatesData.Fs,'yaxis');
    idx_spect=scale_spect(spect);  %calculate index array for spectrogram

    handles.MakeTemplatesData.FreqRows = find((freq >= 1500) & (freq <= 7300));
    %handles.MakeTemplatesData.FreqRows = find((F1 >= 1500) & (F1 <= 7300));
    P1 = idx_spect;
%    P1 = 20*log10(Y1);
    P12 = P1(handles.MakeTemplatesData.FreqRows,:);
    P12 = (P12 - mean(mean(P12)));
    P12 = P12/(sqrt((sum(sum((P12.*P12))))/(size(P12,1) * size(P12,2))));
    Template2 = P12;
    P1 = (P1 - mean(mean(P1)));
    P1 = P1/(sqrt((sum(sum((P1.*P1))))/(size(P1,1) * size(P1,2))));
    Template = P1;
    for j = 2:length(Matches),
        SyllableData2 = handles.MakeTemplatesData.RawData(round(handles.MakeTemplatesData.SyllableBoundaries(Matches(j),1) * handles.MakeTemplatesData.Fs):round(handles.MakeTemplatesData.SyllableBoundaries(Matches(j),2) * handles.MakeTemplatesData.Fs));
        [Y2, F2, T2, P2] = spectrogram(SyllableData2, handles.MakeTemplatesData.FFTWindowLength, handles.MakeTemplatesData.FFTWindowOverlap, handles.MakeTemplatesData.Freq, handles.MakeTemplatesData.Fs, 'yaxis');
        P2 = (P2 - mean(mean(P2)));
        P2 = P2/(sqrt((sum(sum(P2.*P2)))/(size(P2,1) * size(P2,2))));
        [Corr, Lags] = xcorr(sum(P1), sum(P2));
        [X, I] = max(Corr);
        if (Lags(I) > 0)
            if ((size(Template,2) - Lags(I)) <= size(P2,2))
                P2 = P2(:,1:(size(Template,2) - Lags(I)));
                Template(:,((Lags(I)+1):end)) = (Template(:,((Lags(I)+1):end)) + P2)/2;
            else
                Template(:,(Lags(I):(size(P2,2) + Lags(I)))) = (Template(:,(Lags(I):(size(P2,2) + Lags(I)))) + P2)/2;
            end
        else
            if (Lags(I) < 0)
                P2 = P2(:,(abs(Lags(I)) + 1):end);
                if (size(Template,2) <= size(P2,2))
                    P2 = P2(:,1:size(Template,2));
                    Template = (Template + P2)/2;
                else
                    Template(:,(1:size(P2,2))) = (Template(:,(1:size(P2,2))) + P2)/2;
                end
            else
                if (size(Template,2) <= size(P2,2))
                    P2 = P2(:,1:size(Template,2));
                    Template = (Template + P2)/2;
                else
                    Template(:,(1:size(P2,2))) = (Template(:,(1:size(P2,2))) + P2)/2;
                end
            end
        end
        
        [Y2, F2, T2, P2] = spectrogram(SyllableData2, handles.MakeTemplatesData.FFTWindowLength, handles.MakeTemplatesData.FFTWindowOverlap, handles.MakeTemplatesData.Freq, handles.MakeTemplatesData.Fs, 'yaxis');
        P2 = P2(handles.MakeTemplatesData.FreqRows,:);
        P2 = (P2 - mean(mean(P2)));
        P2 = P2/(sqrt((sum(sum(P2.*P2)))/(size(P2,1) * size(P2,2))));
        [Corr, Lags] = xcorr(sum(P12), sum(P2));
        [X, I] = max(Corr);
        if (Lags(I) > 0)
            if ((size(Template2,2) - Lags(I)) <= size(P2,2))
                P2 = P2(:,1:(size(Template2,2) - Lags(I)));
                Template2(:,((Lags(I)+1):end)) = (Template2(:,((Lags(I)+1):end)) + P2)/2;
            else
                Template2(:,(Lags(I):(size(P2,2) + Lags(I)))) = (Template2(:,(Lags(I):(size(P2,2) + Lags(I)))) + P2)/2;
            end
        else
            if (Lags(I) < 0)
                P2 = P2(:,(abs(Lags(I)) + 1):end);
                if (size(Template2,2) <= size(P2,2))
                    P2 = P2(:,1:size(Template2,2));
                    Template2 = (Template2 + P2)/2;
                else
                    Template2(:,(1:size(P2,2))) = (Template2(:,(1:size(P2,2))) + P2)/2;
                end
            else
                if (size(Template2,2) <= size(P2,2))
                    P2 = P2(:,1:size(Template2,2));
                    Template2 = (Template2 + P2)/2;
                else
                    Template2(:,(1:size(P2,2))) = (Template2(:,(1:size(P2,2))) + P2)/2;
                end
            end
        end    
    end
    Templates(i).Template = Template;
    Templates(i).Template2 = Template2;
end 
save(handles.MakeTemplatesData.TemplateFileName, 'Templates');
disp('Saved templates');


function TemplateFileName_Callback(hObject, eventdata, handles)
% hObject    handle to TemplateFileName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TemplateFileName as text
%        str2double(get(hObject,'String')) returns contents of TemplateFileName as a double
handles.MakeTemplatesData.TemplateFileName = get(hObject,'String');
guidata(handles.MakeTemplatesMainWindow, handles);


% --- Executes during object creation, after setting all properties.
function TemplateFileName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TemplateFileName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in QuitButton.
function QuitButton_Callback(hObject, eventdata, handles)
% hObject    handle to QuitButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

FigureTag = findobj('Tag','MakeTemplatesMainWindow');
close(FigureTag);


% --- Executes on button press in ZoomXButton.
function ZoomXButton_Callback(hObject, eventdata, handles)
% hObject    handle to ZoomXButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.ZoomInstructionsText, 'String', 'Left mouse button for left edge, right mouse button for right edge and then press enter (q to quit)');

axes(handles.AmplitudePlot);
TempAxis = axis;
handles.MakeTemplatesData.PresentXCoordinates = TempAxis(1:2);
LeftLine = plot([TempAxis(1) TempAxis(1)], [TempAxis(3) TempAxis(4)], '-.r');
RightLine = plot([TempAxis(2) TempAxis(2)], [TempAxis(3) TempAxis(4)], '-.g');

NewXCoords = TempAxis(1:2);
ZoomBoundaryPatch = patch([NewXCoords(1) NewXCoords(1) NewXCoords(2) NewXCoords(2)], [TempAxis(3) TempAxis(4) TempAxis(4) TempAxis(3)], 'g');
set(ZoomBoundaryPatch, 'FaceAlpha', 0.2);


Flag = 1;

while(Flag)
    [x, y, button] = ginput(1);
    
    if (((length(x) == 0) && (length(y) == 0) && (length(button) == 0)) || (button == 113))
        Flag = 0;
        if (exist('LeftLine', 'var'))
            delete(LeftLine);
        end
        if (exist('RightLine', 'var'))
            delete(RightLine);
        end
        if (exist('ZoomBoundaryPatch', 'var'))
            delete(ZoomBoundaryPatch);
        end
        if (button == 113)
            disp('Quit zoom');
        else
            axes(handles.AmplitudePlot);
            TempAxis(1:2) = handles.MakeTemplatesData.PresentXCoordinates;
            TempAxis(3:4) = handles.MakeTemplatesData.AmpPlotAxis(3:4);
            axis(TempAxis);

            axes(handles.SpectrogramPlot);
            TempAxis(3:4) = handles.MakeTemplatesData.SpecPlotAxis(3:4);
            axis(TempAxis);

            axes(handles.LabelPlot);
            TempAxis(3:4) = handles.MakeTemplatesData.LabelPlotAxis(3:4);
            axis(TempAxis);
        end
            break;
    else
        if (button == 1)
            delete(LeftLine);
            LeftLine = plot([x x], [TempAxis(3) TempAxis(4)], '-.r');
            handles.MakeTemplatesData.PresentXCoordinates(1) = x;
            delete(ZoomBoundaryPatch);
            NewXCoords(1) = x;
            ZoomBoundaryPatch = patch([NewXCoords(1) NewXCoords(1) NewXCoords(2) NewXCoords(2)], [TempAxis(3) TempAxis(4) TempAxis(4) TempAxis(3)], 'g');
            set(ZoomBoundaryPatch, 'FaceAlpha', 0.2);
        else
            if (button == 3)
                delete(RightLine);
                RightLine = plot([x x], [TempAxis(3) TempAxis(4)], '-.g');
                handles.MakeTemplatesData.PresentXCoordinates(2) = x;
                delete(ZoomBoundaryPatch);
                NewXCoords(2) = x;
                ZoomBoundaryPatch = patch([NewXCoords(1) NewXCoords(1) NewXCoords(2) NewXCoords(2)], [TempAxis(3) TempAxis(4) TempAxis(4) TempAxis(3)], 'g');
                set(ZoomBoundaryPatch, 'FaceAlpha', 0.2);                
            end
        end
    end
end

set(handles.ZoomInstructionsText, 'String', 'Watch this space for instructions for zoom and marking syllable boundaries');
guidata(handles.MakeTemplatesMainWindow, handles);


% --- Executes on button press in ZoomYButton.
function ZoomYButton_Callback(hObject, eventdata, handles)
% hObject    handle to ZoomYButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in UnZoomXButton.
function UnZoomXButton_Callback(hObject, eventdata, handles)
% hObject    handle to UnZoomXButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.MakeTemplatesData.PresentXCoordinates = handles.MakeTemplatesData.AmpPlotAxis(1:2);
guidata(handles.MakeTemplatesMainWindow, handles);

axes(handles.AmplitudePlot);
axis(handles.MakeTemplatesData.AmpPlotAxis);

axes(handles.SpectrogramPlot);
axis(handles.MakeTemplatesData.SpecPlotAxis);

axes(handles.LabelPlot);
axis(handles.MakeTemplatesData.LabelPlotAxis);

% --- Executes on button press in UnZoomYButton.
function UnZoomYButton_Callback(hObject, eventdata, handles)
% hObject    handle to UnZoomYButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

axes(handles.AmplitudePlot);
axis(handles.MakeTemplatesData.AmpPlotAxis);

axes(handles.SpectrogramPlot);
axis(handles.MakeTemplatesData.SpecPlotAxis);

axes(handles.LabelPlot);
axis(handles.MakeTemplatesData.LabelPlotAxis);
