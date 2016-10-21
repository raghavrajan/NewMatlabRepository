function varargout = ASSLAdjustSyllableBoundaries(varargin)
% ASSLADJUSTSYLLABLEBOUNDARIES MATLAB code for ASSLAdjustSyllableBoundaries.fig
%      ASSLADJUSTSYLLABLEBOUNDARIES, by itself, creates a new ASSLADJUSTSYLLABLEBOUNDARIES or raises the existing
%      singleton*.
%
%      H = ASSLADJUSTSYLLABLEBOUNDARIES returns the handle to a new ASSLADJUSTSYLLABLEBOUNDARIES or the handle to
%      the existing singleton*.
%
%      ASSLADJUSTSYLLABLEBOUNDARIES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ASSLADJUSTSYLLABLEBOUNDARIES.M with the given input arguments.
%
%      ASSLADJUSTSYLLABLEBOUNDARIES('Property','Value',...) creates a new ASSLADJUSTSYLLABLEBOUNDARIES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ASSLAdjustSyllableBoundaries_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ASSLAdjustSyllableBoundaries_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ASSLAdjustSyllableBoundaries

% Last Modified by GUIDE v2.5 22-Feb-2013 03:05:23

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ASSLAdjustSyllableBoundaries_OpeningFcn, ...
                   'gui_OutputFcn',  @ASSLAdjustSyllableBoundaries_OutputFcn, ...
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


% --- Executes just before ASSLAdjustSyllableBoundaries is made visible.
function ASSLAdjustSyllableBoundaries_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ASSLAdjustSyllableBoundaries (see VARARGIN)

% Choose default command line output for ASSLAdjustSyllableBoundaries
handles.output = hObject;

% Initialise some variables using either default values or values that have
% been passed while calling the function (RR 07 Jan 2013)

if (nargin >= 1)
    Temp = varargin{1};
    handles.ASB = Temp.ASSL;
    TempLabels = [];
    handles.ASB.AllSyllIndices = [];
    handles.ASB.AllSyllLabels = [];
    if (isfield(handles.ASB, 'SyllLabels'))
        for i = 1:length(handles.ASB.SyllLabels),
            TempLabels = [TempLabels handles.ASB.SyllLabels{i}];
            handles.ASB.AllSyllIndices = [handles.ASB.AllSyllIndices; [ones(size(handles.ASB.SyllLabels{i}'))*i (1:1:length(handles.ASB.SyllLabels{i}))']];
            handles.ASB.AllSyllLabels = [handles.ASB.AllSyllLabels; [handles.ASB.SyllLabels{i}']];
        end
        TempLabels = unique(TempLabels);
        %Temp = double(TempLabels);
        %TempLabels = TempLabels(find((Temp >=97) & (Temp <= 122)));
        for i = 1:length(TempLabels),
            handles.ASB.UniqueSyllLabels{i} = TempLabels(i);
        end
        handles.ASB.SyllCounter = 1;
        set(handles.SyllableListBox, 'String', handles.ASB.UniqueSyllLabels);
        set(handles.SyllableListBox, 'Value', handles.ASB.SyllCounter);
        
        handles.ASB.SyllExampleCounter = 1;
        
        handles.ASB.SyllData = ASSLGetIndividualSyllData(handles);
        
        % Variables for plotting
%         Dir = handles.ASB.DirName;
%         FileName = handles.ASB.DirName{handles.ASB.SyllCounter}.FileName{handles.SSASA.SyllExampleCounter};
% 
%         FileType = handles.SSASA.FileType;
%         SpecAxis = handles.ASSLExampleSpectAxes;
%         SpikeTrainAxis = handles.ASSLExampleAmplitudeAxes;
%         SpikeTimes = handles.SSASA.SyllData{handles.SSASA.SyllCounter}.SpikeTimes{handles.SSASA.SyllExampleCounter}; 
%         TimeRange = handles.SSASA.SyllData{handles.SSASA.SyllCounter}.SyllDur{handles.SSASA.SyllExampleCounter};
%         TimeRange(1) = TimeRange(1) - 0.2;
%         TimeRange(2) = TimeRange(2) + 0.2;
%         FileDur = handles.SSASA.SyllData{handles.SSASA.SyllCounter}.FileDur{handles.SSASA.SyllExampleCounter};
%         
%         SSAPlotExampleSpectST(Dir, FileName, FileType, SpecAxis, SpikeTrainAxis, SpikeTimes, TimeRange, FileDur, 0.2, 0.2);
%         
        axes(handles.ASSLAllSyllWaveformAxes);
        cla(handles.ASSLAllSyllWaveformAxes);
        ASSLPlotSyllAmplitudeWaveforms(handles.ASB.SyllData{handles.ASB.SyllCounter}, 'NormalOnset');
        
        axes(handles.ASSLAllSyllOffsetWaveformAxes);
        cla(handles.ASSLAllSyllOffsetWaveformAxes);
        ASSLPlotSyllAmplitudeWaveforms(handles.ASB.SyllData{handles.ASB.SyllCounter}, 'NormalOffset');
    else
        set(handles.SyllableListBox, 'String', []);
    end
    
end
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ASSLAdjustSyllableBoundaries wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = ASSLAdjustSyllableBoundaries_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in SyllableListBox.
function SyllableListBox_Callback(hObject, eventdata, handles)
% hObject    handle to SyllableListBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns SyllableListBox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from SyllableListBox
% Variables for plotting
handles.ASB.SyllCounter = get(hObject, 'Value');
handles.ASB.SyllExampleCounter = 1;

% Dir = handles.SSASA.RawDataDirectory;
% FileName = handles.SSASA.SyllData{handles.SSASA.SyllCounter}.FileName{handles.SSASA.SyllExampleCounter};
% FileType = handles.SSASA.FileType;
% SpecAxis = handles.ASSLExampleSpectAxes;
% SpikeTrainAxis = handles.ASSLExampleAmplitudeAxes;
% SpikeTimes = handles.SSASA.SyllData{handles.SSASA.SyllCounter}.SpikeTimes{handles.SSASA.SyllExampleCounter}; 
% TimeRange = handles.SSASA.SyllData{handles.SSASA.SyllCounter}.SyllDur{handles.SSASA.SyllExampleCounter};
% TimeRange(1) = TimeRange(1) - 0.2;
% TimeRange(2) = TimeRange(2) + 0.2;
% FileDur = handles.SSASA.SyllData{handles.SSASA.SyllCounter}.FileDur{handles.SSASA.SyllExampleCounter};
% 
% SSAPlotExampleSpectST(Dir, FileName, FileType, SpecAxis, SpikeTrainAxis, SpikeTimes, TimeRange, FileDur, 0.2, 0.2);


axes(handles.ASSLAllSyllWaveformAxes);
cla(handles.ASSLAllSyllWaveformAxes);
ASSLPlotSyllAmplitudeWaveforms(handles.ASB.SyllData{handles.ASB.SyllCounter}, 'NormalOnset');

axes(handles.ASSLAllSyllOffsetWaveformAxes);
cla(handles.ASSLAllSyllOffsetWaveformAxes);
ASSLPlotSyllAmplitudeWaveforms(handles.ASB.SyllData{handles.ASB.SyllCounter}, 'NormalOffset');

if (isfield(handles.ASB.SyllData{1}, 'Adjust'))
    axes(handles.ASSLAdjustedSyllWaveformAxes);
    cla(handles.ASSLAdjustedSyllWaveformAxes);
    ASSLPlotSyllAmplitudeWaveforms(handles.ASB.SyllData{handles.ASB.SyllCounter}, 'AdjustedOnset');

    axes(handles.ASSLAdjustedSyllOffsetWaveformAxes);
    cla(handles.ASSLAdjustedSyllOffsetWaveformAxes);
    ASSLPlotSyllAmplitudeWaveforms(handles.ASB.SyllData{handles.ASB.SyllCounter}, 'AdjustedOffset');
end

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function SyllableListBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SyllableListBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
-15

% --- Executes on button press in AdjustWithCorrelationsButton.
function AdjustWithCorrelationsButton_Callback(hObject, eventdata, handles)
% hObject    handle to AdjustWithCorrelationsButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

for j = 1:length(handles.ASB.SyllData),
    SyllableStart1 = find(handles.ASB.SyllData{j}.Time{1} >= 0, 1, 'first');
    SyllableEnd1 = find(handles.ASB.SyllData{j}.Time{1} >= handles.ASB.SyllData{j}.OffsetTime{1}, 1, 'first');
    for i = 2:length(handles.ASB.SyllData{j}.Time),
        SyllableStart = find(handles.ASB.SyllData{j}.Time{i} >= 0, 1, 'first');
        [c, lags] = xcorr(handles.ASB.SyllData{j}.Amplitude{1}(1:SyllableStart1+160), handles.ASB.SyllData{j}.Amplitude{i}(1:SyllableStart+160));
        [maxval, maxind] = max(c);
        handles.ASB.SyllData{j}.Adjust(i,1) = lags(maxind)/32;

        handles.ASB.SyllData{j}.Adjust(i,1) = lags(maxind)/32;
        handles.ASB.AdjSyllOnsets{handles.ASB.SyllData{j}.FileSyllIndices(i,1)}(handles.ASB.SyllData{j}.FileSyllIndices(i,2)) = handles.ASB.SyllOnsets{handles.ASB.SyllData{j}.FileSyllIndices(i,1)}(handles.ASB.SyllData{j}.FileSyllIndices(i,2)) + handles.ASB.SyllData{j}.Adjust(i,1);
        
        SyllableEnd = find(handles.ASB.SyllData{j}.Time{i} >= handles.ASB.SyllData{j}.OffsetTime{i}, 1, 'first');
        [c, lags] = xcorr(handles.ASB.SyllData{j}.Amplitude{1}(SyllableEnd1-160:end), handles.ASB.SyllData{j}.Amplitude{i}(SyllableEnd-160:end));
        [maxval, maxind] = max(c);
        handles.ASB.SyllData{j}.Adjust(i,2) = (lags(maxind))/32;
        
        handles.ASB.AdjSyllOffsets{handles.ASB.SyllData{j}.FileSyllIndices(i,1)}(handles.ASB.SyllData{j}.FileSyllIndices(i,2)) = handles.ASB.SyllOffsets{handles.ASB.SyllData{j}.FileSyllIndices(i,1)}(handles.ASB.SyllData{j}.FileSyllIndices(i,2)) + handles.ASB.SyllData{j}.Adjust(i,2);
    end
end

guidata(hObject, handles);

% --- Executes on button press in ThresholdAdjustButton.
function ThresholdAdjustButton_Callback(hObject, eventdata, handles)
% hObject    handle to ThresholdAdjustButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Threshold = inputdlg('Set new threshold', 'Threshold box');
Threshold = str2double(Threshold{1});

for j = 1:length(handles.ASB.SyllData),
    for i = 1:length(handles.ASB.SyllData{j}.Time),
        SyllableStart = find(handles.ASB.SyllData{j}.Time{i} >= 0, 1, 'first');
        ThresholdCrossing = find(handles.ASB.SyllData{j}.Amplitude{i}(1:SyllableStart+160) <= Threshold, 1, 'last');
        if (isempty(ThresholdCrossing))
            handles.ASB.SyllData{j}.Adjust(i,1) = 0;
        else
            handles.ASB.SyllData{j}.Adjust(i,1) = (ThresholdCrossing - SyllableStart)/32;
        end
        
        handles.ASB.AdjSyllOnsets{handles.ASB.SyllData{j}.FileSyllIndices(i,1)}(handles.ASB.SyllData{j}.FileSyllIndices(i,2)) = handles.ASB.SyllOnsets{handles.ASB.SyllData{j}.FileSyllIndices(i,1)}(handles.ASB.SyllData{j}.FileSyllIndices(i,2)) + handles.ASB.SyllData{j}.Adjust(i,1);
        
        SyllableEnd = find(handles.ASB.SyllData{j}.Time{i} >= handles.ASB.SyllData{j}.OffsetTime{i}, 1, 'first');
        ThresholdCrossing = find(handles.ASB.SyllData{j}.Amplitude{i}(SyllableEnd-160:end) <= Threshold, 1, 'first');
        ThresholdCrossing = ThresholdCrossing + SyllableEnd - 160 - 1;
        
        if (isempty(ThresholdCrossing))
            handles.ASB.SyllData{j}.Adjust(i,2) = 0;
        else
            handles.ASB.SyllData{j}.Adjust(i,2) = (ThresholdCrossing - SyllableEnd)/32;
        end
        
        handles.ASB.AdjSyllOffsets{handles.ASB.SyllData{j}.FileSyllIndices(i,1)}(handles.ASB.SyllData{j}.FileSyllIndices(i,2)) = handles.ASB.SyllOffsets{handles.ASB.SyllData{j}.FileSyllIndices(i,1)}(handles.ASB.SyllData{j}.FileSyllIndices(i,2)) + handles.ASB.SyllData{j}.Adjust(i,2);
    end
end

guidata(hObject, handles);

% --- Executes on button press in FirstDerivAdjustButton.
function FirstDerivAdjustButton_Callback(hObject, eventdata, handles)
% hObject    handle to FirstDerivAdjustButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

SmoothingWindow = ones(1,round(2*32))/round(2*32);

for j = 1:length(handles.ASB.SyllData),
    for i = 1:length(handles.ASB.SyllData{j}.Time),
        FirstDerivative = diff(handles.ASB.SyllData{j}.Amplitude{i});
        FirstDerivative = conv(FirstDerivative, SmoothingWindow, 'same');
        
        SyllableStart = find(handles.ASB.SyllData{j}.Time{i} >= 0, 1, 'first');
        [MaxVal, ThresholdCrossing] = max(FirstDerivative(1:SyllableStart+160));
        if (isempty(ThresholdCrossing))
            handles.ASB.SyllData{j}.Adjust(i,1) = 0;
        else
            handles.ASB.SyllData{j}.Adjust(i,1) = (ThresholdCrossing - SyllableStart)/32;
        end
        
        handles.ASB.AdjSyllOnsets{handles.ASB.SyllData{j}.FileSyllIndices(i,1)}(handles.ASB.SyllData{j}.FileSyllIndices(i,2)) = handles.ASB.SyllOnsets{handles.ASB.SyllData{j}.FileSyllIndices(i,1)}(handles.ASB.SyllData{j}.FileSyllIndices(i,2)) + handles.ASB.SyllData{j}.Adjust(i,1);
        
        SyllableEnd = find(handles.ASB.SyllData{j}.Time{i} >= handles.ASB.SyllData{j}.OffsetTime{i}, 1, 'first');
        [MaxVal, ThresholdCrossing] = max(-FirstDerivative(SyllableEnd-160:end));
        ThresholdCrossing = ThresholdCrossing + SyllableEnd - 160 - 1;
        
        if (isempty(ThresholdCrossing))
            handles.ASB.SyllData{j}.Adjust(i,2) = 0;
        else
            handles.ASB.SyllData{j}.Adjust(i,2) = (ThresholdCrossing - SyllableEnd)/32;
        end
        
        handles.ASB.AdjSyllOffsets{handles.ASB.SyllData{j}.FileSyllIndices(i,1)}(handles.ASB.SyllData{j}.FileSyllIndices(i,2)) = handles.ASB.SyllOffsets{handles.ASB.SyllData{j}.FileSyllIndices(i,1)}(handles.ASB.SyllData{j}.FileSyllIndices(i,2)) + handles.ASB.SyllData{j}.Adjust(i,2);
    end
end

disp('Finished adjusting boundaries with first derivative peaks');
guidata(hObject, handles);

% --- Executes on button press in SecondDerivAdjustButton.
function SecondDerivAdjustButton_Callback(hObject, eventdata, handles)
% hObject    handle to SecondDerivAdjustButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

SmoothingWindow = ones(1,round(5*32))/round(5*32);

for j = 1:length(handles.ASB.SyllData),
    for i = 1:length(handles.ASB.SyllData{j}.Time),
        SecondDerivative = diff(handles.ASB.SyllData{j}.Amplitude{i}, 2);
        SecondDerivative = conv(SecondDerivative, SmoothingWindow, 'same');
        
        SyllableStart = find(handles.ASB.SyllData{j}.Time{i} >= 0, 1, 'first');
        [MaxVal, ThresholdCrossing] = max(SecondDerivative(1:SyllableStart+160));
        if (isempty(ThresholdCrossing))
            handles.ASB.SyllData{j}.Adjust(i,1) = 0;
        else
            handles.ASB.SyllData{j}.Adjust(i,1) = (ThresholdCrossing - SyllableStart)/32;
        end
        
        handles.ASB.AdjSyllOnsets{handles.ASB.SyllData{j}.FileSyllIndices(i,1)}(handles.ASB.SyllData{j}.FileSyllIndices(i,2)) = handles.ASB.SyllOnsets{handles.ASB.SyllData{j}.FileSyllIndices(i,1)}(handles.ASB.SyllData{j}.FileSyllIndices(i,2)) + handles.ASB.SyllData{j}.Adjust(i,1);
        
        SyllableEnd = find(handles.ASB.SyllData{j}.Time{i} >= handles.ASB.SyllData{j}.OffsetTime{i}, 1, 'first');
        [MaxVal, ThresholdCrossing] = max(-SecondDerivative(SyllableEnd-160:end));
        ThresholdCrossing = ThresholdCrossing + SyllableEnd - 160 - 1;
        
        if (isempty(ThresholdCrossing))
            handles.ASB.SyllData{j}.Adjust(i,2) = 0;
        else
            handles.ASB.SyllData{j}.Adjust(i,2) = (ThresholdCrossing - SyllableEnd)/32;
        end
        
        handles.ASB.AdjSyllOffsets{handles.ASB.SyllData{j}.FileSyllIndices(i,1)}(handles.ASB.SyllData{j}.FileSyllIndices(i,2)) = handles.ASB.SyllOffsets{handles.ASB.SyllData{j}.FileSyllIndices(i,1)}(handles.ASB.SyllData{j}.FileSyllIndices(i,2)) + handles.ASB.SyllData{j}.Adjust(i,2);
    end
end

disp('Finished adjusting boundaries with second derivative peaks');
guidata(hObject, handles);

% --- Executes on button press in FWHMAdjustButton.
function FWHMAdjustButton_Callback(hObject, eventdata, handles)
% hObject    handle to FWHMAdjustButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

for j = 1:length(handles.ASB.SyllData),
    for i = 1:length(handles.ASB.SyllData{j}.Time),
        Threshold = (max(handles.ASB.SyllData{j}.Amplitude{i}) + min(handles.ASB.SyllData{j}.Amplitude{i}))/2;
        
        ThresholdCrossings = FindThresholdCrossings(handles.ASB.SyllData{j}.Amplitude{i}, Threshold);
        
        SyllableStart = find(handles.ASB.SyllData{j}.Time{i} >= 0, 1, 'first');
        SyllableEnd = find(handles.ASB.SyllData{j}.Time{i} >= handles.ASB.SyllData{j}.OffsetTime{i}, 1, 'first');
        
        if (isempty(ThresholdCrossings))
            handles.ASB.SyllData{j}.Adjust(i,1) = 0;
            handles.ASB.SyllData{j}.Adjust(i,2) = 0;
        else
            if (handles.ASB.SyllData{j}.Amplitude{i}(SyllableStart) < Threshold)
                MinIndex = (SyllableStart - 1) + find(handles.ASB.SyllData{j}.Amplitude{i}(SyllableStart:end) >= Threshold, 1, 'first');
            else
                MinIndex = find(handles.ASB.SyllData{j}.Amplitude{i}(1:SyllableStart) <= Threshold, 1, 'last');
            end
            
            if (isempty(MinIndex))
                handles.ASB.SyllData{j}.Adjust(i,1) = 0;
            else
                handles.ASB.SyllData{j}.Adjust(i,1) = handles.ASB.SyllData{j}.Time{i}(MinIndex);
            end
            
            if (handles.ASB.SyllData{j}.Amplitude{i}(SyllableEnd) > Threshold)
                MinIndex = (SyllableEnd - 1) + find(handles.ASB.SyllData{j}.Amplitude{i}(SyllableEnd:end) <= Threshold, 1, 'first');
            else
                MinIndex = find(handles.ASB.SyllData{j}.Amplitude{i}(1:SyllableEnd) >= Threshold, 1, 'last');
            end
            
            if (isempty(MinIndex))
                handles.ASB.SyllData{j}.Adjust(i,2) = 0;
            else
                handles.ASB.SyllData{j}.Adjust(i,2) = handles.ASB.SyllData{j}.Time{i}(MinIndex) - handles.ASB.SyllData{j}.OffsetTime{i};
            end
        end
        
        handles.ASB.AdjSyllOnsets{handles.ASB.SyllData{j}.FileSyllIndices(i,1)}(handles.ASB.SyllData{j}.FileSyllIndices(i,2)) = handles.ASB.SyllOnsets{handles.ASB.SyllData{j}.FileSyllIndices(i,1)}(handles.ASB.SyllData{j}.FileSyllIndices(i,2)) + handles.ASB.SyllData{j}.Adjust(i,1);
        handles.ASB.AdjSyllOffsets{handles.ASB.SyllData{j}.FileSyllIndices(i,1)}(handles.ASB.SyllData{j}.FileSyllIndices(i,2)) = handles.ASB.SyllOffsets{handles.ASB.SyllData{j}.FileSyllIndices(i,1)}(handles.ASB.SyllData{j}.FileSyllIndices(i,2)) + handles.ASB.SyllData{j}.Adjust(i,2);
    end
end

guidata(hObject, handles);

disp('Finished adjusting using full-width-at-half-maximum as threshold');

% --- Executes on button press in ReturnAdjBoundariesButton.
function ReturnAdjBoundariesButton_Callback(hObject, eventdata, handles)
% hObject    handle to ReturnAdjBoundariesButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ASSLMainWindow = findobj('Name', 'AutoSongSegmentLabel');
if (length(ASSLMainWindow) > 0)
    ASSLData = guidata(ASSLMainWindow);
    ASSLData.ASSL.AdjSyllOnsets = handles.ASB.AdjSyllOnsets;
    ASSLData.ASSL.AdjSyllOffsets = handles.ASB.AdjSyllOffsets;
    guidata(ASSLMainWindow, ASSLData);
    msgbox('Returned adjusted syllable boundaries to Auto Song Segment Label');
end
guidata(hObject, handles);
