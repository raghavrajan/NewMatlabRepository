function [] = ASSLTMPFPlotData(handles)

cla(handles.ASSLTemplateMatchFeaturePlotAxis);
axes(handles.ASSLTemplateMatchFeaturePlotAxis);
plot(handles.ASSLTMPF.FeatVals(:,handles.ASSLTMPF.XVal), handles.ASSLTMPF.FeatVals(:,handles.ASSLTMPF.YVal), 'k+', 'MarkerSize', 2);
hold on;
plot(handles.ASSLTMPF.FeatVals(handles.ASSLTMPF.SyllIndex, handles.ASSLTMPF.XVal), handles.ASSLTMPF.FeatVals(handles.ASSLTMPF.SyllIndex, handles.ASSLTMPF.YVal), 'ro', 'MarkerSize', 6);
xlabel(handles.ASSLTMPF.FeatNames{handles.ASSLTMPF.XVal}, 'FontSize', 12, 'FontWeight', 'bold');
ylabel(handles.ASSLTMPF.FeatNames{handles.ASSLTMPF.YVal}, 'FontSize', 12, 'FontWeight', 'bold');
set(gca, 'FontSize', 12, 'FontWeight', 'bold');
axis tight;

SyllLabels = get(handles.TemplateListBox, 'String');
Indices = find(handles.DataStruct.SyllIndexLabels == SyllLabels{handles.ASSLTMPF.TemplateIndex});
plot(handles.ASSLTMPF.FeatVals(Indices, handles.ASSLTMPF.XVal), handles.ASSLTMPF.FeatVals(Indices, handles.ASSLTMPF.YVal), 'r+', 'MarkerSize', 2);

cla(handles.ASSLTMPFSpecAxis);
set(handles.SyllableNoText, 'String', ['Syllable # : ', num2str(handles.ASSLTMPF.SyllIndex), ' of ', num2str(size(handles.ASSLTMPF.FeatVals,1)), ' syllables']);

[RawData, Fs] = ASSLGetRawData(handles.DataStruct.DirName, handles.DataStruct.FileName{handles.DataStruct.SyllIndices(handles.ASSLTMPF.SyllIndex,1)}, handles.DataStruct.FileType, handles.DataStruct.SongChanNo);
    
Time = (1:1:length(RawData))/Fs;
[LogAmplitude] = ASSLCalculateLogAmplitude(RawData, Fs, Time, handles.DataStruct.FFTWinSizeSegmenting, handles.DataStruct.FFTWinOverlapSegmenting);

TimeRange(1) = handles.DataStruct.SyllOnsets{handles.DataStruct.SyllIndices(handles.ASSLTMPF.SyllIndex,1)}(handles.DataStruct.SyllIndices(handles.ASSLTMPF.SyllIndex,2))/1000 - 0.005;
TimeRange(2) = handles.DataStruct.SyllOffsets{handles.DataStruct.SyllIndices(handles.ASSLTMPF.SyllIndex,1)}(handles.DataStruct.SyllIndices(handles.ASSLTMPF.SyllIndex,2))/1000 + 0.005;

ASSLPlotData(handles.DataStruct.DirName, handles.DataStruct.FileName{handles.DataStruct.SyllIndices(handles.ASSLTMPF.SyllIndex,1)}, handles.DataStruct.FileType, Time, LogAmplitude, handles.ASSLTMPFSpecAxis, handles.ASSLTMPFAmpAxis, handles.DataStruct.Threshold{handles.DataStruct.SyllIndices(handles.ASSLTMPF.SyllIndex,1)}, (TimeRange(1)+0.005)*1000, (TimeRange(2)-0.005)*1000, [], TimeRange);

cla(handles.TemplateSpecAxis);
axes(handles.TemplateSpecAxis);
contourf(handles.DataStruct.Templates{1}.SyllableTemplates{handles.ASSLTMPF.TemplateIndex}{1}.MotifTemplate(4).MotifTemplate);
colormap('jet');

cla(handles.TemplateMatchHistAxis);
axes(handles.TemplateMatchHistAxis);
%Edges = 0:0.1:max(max(handles.DataStruct.TemplateMatchValues));
%plot(Edges, histc(handles.DataStruct.TemplateMatchValues(:,handles.ASSLTMPF.TemplateIndex), Edges), 'k');
%axis tight;