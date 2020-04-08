function [] = ASSLPlotSyllFFBoundaries(DataStruct)

FF = DataStruct.DataStruct.FeatValues(:,end);

if (isfield(DataStruct.ASSLCSFFB, 'FFSyllBoundaries'))
    FFSyllBoundaries = DataStruct.ASSLCSFFB.FFSyllBoundaries;
else
    FFSyllBoundaries = [];
end

Indices = find(DataStruct.DataStruct.SyllIndexLabels == DataStruct.ASSLCSFFB.UniqueSyllLabels(DataStruct.ASSLCSFFB.SyllIndex));

ElapsedTime = 0;

cla(DataStruct.axes1);
for i = 1:length(Indices),
    if (i < DataStruct.ASSLCSFFB.CurrentStartIndex)
        continue;
    end
    
    if (i >= (DataStruct.ASSLCSFFB.CurrentStartIndex + DataStruct.ASSLCSFFB.NumExamples))
        continue;
    end
    
    [Syllable, Fs] = ASSLGetRawData(DataStruct.DataStruct.FileDir{DataStruct.DataStruct.SyllIndices(Indices(i), 1)}, DataStruct.DataStruct.FileName{DataStruct.DataStruct.SyllIndices(Indices(i), 1)}, DataStruct.DataStruct.FileType, DataStruct.DataStruct.SongChanNo);
    if (isfield(DataStruct.DataStruct, 'AdjSyllOnsets'))
        StartTime = (ceil(Fs * DataStruct.DataStruct.AdjSyllOnsets{DataStruct.DataStruct.SyllIndices(Indices(i), 1)}(DataStruct.DataStruct.SyllIndices(Indices(i),2))/1000));
    else
        StartTime = (ceil(Fs * DataStruct.DataStruct.SyllOnsets{DataStruct.DataStruct.SyllIndices(Indices(i), 1)}(DataStruct.DataStruct.SyllIndices(Indices(i),2))/1000));
    end
    if (StartTime == 0)
        StartTime = 1;
    else
        if (StartTime >= length(Syllable))
            StartTime = (length(Syllable) - 1);
        end
    end
    
    if (isfield(DataStruct.DataStruct, 'AdjSyllOnsets'))
        EndTime = (ceil(Fs * DataStruct.DataStruct.AdjSyllOffsets{DataStruct.DataStruct.SyllIndices(Indices(i), 1)}(DataStruct.DataStruct.SyllIndices(Indices(i),2))/1000));
    else
        EndTime = (ceil(Fs * DataStruct.DataStruct.SyllOffsets{DataStruct.DataStruct.SyllIndices(Indices(i), 1)}(DataStruct.DataStruct.SyllIndices(Indices(i),2))/1000));
    end
    if (EndTime == 0)
        EndTime = 2;
    else
        if (EndTime >= length(Syllable))
            EndTime = length(Syllable);
        end
    end
    
    Syllable = Syllable(StartTime:EndTime);
    Time = ((1:1:length(Syllable))/Fs) + (ElapsedTime) + 0.04;
    ElapsedTime = Time(end);
    
    TempFF = DataStruct.DataStruct.Raw.FundamentalFrequency{DataStruct.DataStruct.SyllIndices(Indices(i), 1)}{DataStruct.DataStruct.SyllIndices(Indices(i),2)};
    TempFF_x = linspace(Time(1), Time(end), length(TempFF));
    StartIndex = find(TempFF_x <= (Time(1) + FFSyllBoundaries(Indices(i), 1)/Fs), 1, 'last');
    if (isempty(StartIndex))
        StartIndex = 1;
    end

    EndIndex = find(TempFF_x >= (Time(1) + FFSyllBoundaries(Indices(i), 2)/Fs), 1, 'first');
    if (isempty(EndIndex))
        EndIndex = length(TempFF_x);
    end

    if (~isempty(TempFF))
        FF(Indices(i)) = mean(TempFF(StartIndex:EndIndex));
    else
        FF(Indices(i)) = NaN;
    end

    PlotSpectrogramInAxis_SongVar(Syllable, Time, Fs, DataStruct.axes1)
    axes(DataStruct.axes1);
    hold on;
    axis([0 1.1*Time(end) 300 8000]);

    plot(TempFF_x, TempFF, 'b', 'LineWidth', 2); 
    LeftBoundaryLine{i} = plot([(Time(1) + FFSyllBoundaries(Indices(i),1)/Fs) (Time(1) + FFSyllBoundaries(Indices(i),1)/Fs)], [300 8000], 'r--', 'LineWidth', 2);
    RightBoundaryLine{i} = plot([(Time(1) + FFSyllBoundaries(Indices(i),2)/Fs) (Time(1) + FFSyllBoundaries(Indices(i),2)/Fs)], [300 8000], 'g--', 'LineWidth', 2);
    plot(TempFF_x(StartIndex:EndIndex), FF(Indices(i)), 'g--', 'LineWidth', 2);
end

set(DataStruct.SyllableIdentityLabel, 'String', ['Syllable ', DataStruct.ASSLCSFFB.UniqueSyllLabels(DataStruct.ASSLCSFFB.SyllIndex), ': #', num2str(DataStruct.ASSLCSFFB.SyllIndex), ' of ', num2str(length(DataStruct.ASSLCSFFB.UniqueSyllLabels)), ' syllables; Time dur of FF = ', num2str(TempFF_x(EndIndex) - TempFF_x(StartIndex)), ' sec; ', num2str(DataStruct.ASSLCSFFB.CurrentStartIndex), ':', num2str(DataStruct.ASSLCSFFB.CurrentStartIndex + DataStruct.ASSLCSFFB.NumExamples - 1), ' of ', num2str(length(Indices))]);

disp('Finished plotting boundaries');
