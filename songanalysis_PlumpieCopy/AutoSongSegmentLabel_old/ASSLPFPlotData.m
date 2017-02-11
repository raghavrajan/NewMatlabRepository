function [] = ASSLPFPlotData(handles)

cla(handles.ASSLFeaturePlotAxis);
axes(handles.ASSLFeaturePlotAxis);
plot(handles.ASSLPF.FeatVals(:,handles.ASSLPF.XVal), handles.ASSLPF.FeatVals(:,handles.ASSLPF.YVal), 'k+', 'MarkerSize', 2);
hold on;
plot(handles.ASSLPF.FeatVals(handles.ASSLPF.SyllIndex, handles.ASSLPF.XVal), handles.ASSLPF.FeatVals(handles.ASSLPF.SyllIndex, handles.ASSLPF.YVal), 'ro', 'MarkerSize', 6);
xlabel(handles.ASSLPF.FeatNames{handles.ASSLPF.XVal}, 'FontSize', 12, 'FontWeight', 'bold');
ylabel(handles.ASSLPF.FeatNames{handles.ASSLPF.YVal}, 'FontSize', 12, 'FontWeight', 'bold');
set(gca, 'FontSize', 12, 'FontWeight', 'bold');
axis tight;

cla(handles.ASSLPFSpecAxis);
set(handles.SyllableNoText, 'String', ['Syllable # : ', num2str(handles.ASSLPF.SyllIndex), ' of ', num2str(size(handles.ASSLPF.FeatVals,1)), ' syllables']);

[RawData, Fs] = ASSLGetRawData(handles.DataStruct.DirName, handles.DataStruct.FileName{handles.DataStruct.SyllIndices(handles.ASSLPF.SyllIndex,1)}, handles.DataStruct.FileType, handles.DataStruct.SongChanNo);
    
Time = (1:1:length(RawData))/Fs;
[LogAmplitude] = ASSLCalculateLogAmplitude(RawData, Fs, Time, handles.DataStruct.FFTWinSizeSegmenting, handles.DataStruct.FFTWinOverlapSegmenting);

TimeRange(1) = handles.DataStruct.SyllOnsets{handles.DataStruct.SyllIndices(handles.ASSLPF.SyllIndex,1)}(handles.DataStruct.SyllIndices(handles.ASSLPF.SyllIndex,2))/1000 - 0.005;
TimeRange(2) = handles.DataStruct.SyllOffsets{handles.DataStruct.SyllIndices(handles.ASSLPF.SyllIndex,1)}(handles.DataStruct.SyllIndices(handles.ASSLPF.SyllIndex,2))/1000 + 0.005;

ASSLPlotData(handles.DataStruct.DirName, handles.DataStruct.FileName{handles.DataStruct.SyllIndices(handles.ASSLPF.SyllIndex,1)}, handles.DataStruct.FileType, Time, LogAmplitude, handles.ASSLPFSpecAxis, handles.ASSLPFAmpAxis, handles.DataStruct.Threshold{handles.DataStruct.SyllIndices(handles.ASSLPF.SyllIndex,1)}, (TimeRange(1)+0.005)*1000, (TimeRange(2)-0.005)*1000, [], TimeRange);
