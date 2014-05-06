function [FFSyllBoundaries, FF] = CalculatePlotSyllFFBoundaries(DataStruct)

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
    [Syllable, Fs] = ASSLGetRawData(DataStruct.DataStruct.DirName, DataStruct.DataStruct.FileName{DataStruct.DataStruct.SyllIndices(Indices(i), 1)}, DataStruct.DataStruct.FileType, DataStruct.DataStruct.SongChanNo);
    StartTime = (ceil(Fs * DataStruct.DataStruct.SyllOnsets{DataStruct.DataStruct.SyllIndices(Indices(i), 1)}(DataStruct.DataStruct.SyllIndices(Indices(i),2))/1000));
    if (StartTime == 0)
        StartTime = 1;
    else
        if (StartTime >= length(Syllable))
            StartTime = (length(Syllable) - 1);
        end
    end
    EndTime = (ceil(Fs * DataStruct.DataStruct.SyllOffsets{DataStruct.DataStruct.SyllIndices(Indices(i), 1)}(DataStruct.DataStruct.SyllIndices(Indices(i),2))/1000));
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
    
    if (strfind(DataStruct.ASSLCSFFB.BoundaryChoices{DataStruct.ASSLCSFFB.BoundaryChoice}, 'Specify limits by percent from beginning'))
        % Two values; left and right boundaries considering 0 as the start time of
        % the syllable. Boundaries are calculated and stored as indices
        
        FFSyllBoundaries(Indices(i), 1) = ceil(DataStruct.ASSLCSFFB.StartLimit/100 * length(Syllable));
        FFSyllBoundaries(Indices(i), 2) = ceil(DataStruct.ASSLCSFFB.EndLimit/100 * length(Syllable));
    else
        if (strfind(DataStruct.ASSLCSFFB.BoundaryChoices{DataStruct.ASSLCSFFB.BoundaryChoice}, 'Specify limits by percent from end'))
            FFSyllBoundaries(Indices(i), 1) = length(Syllable) - ceil(DataStruct.ASSLCSFFB.StartLimit/100 * length(Syllable));
            FFSyllBoundaries(Indices(i), 2) = length(Syllable) - ceil(DataStruct.ASSLCSFFB.EndLimit/100 * length(Syllable));
        else
            if (strfind(DataStruct.ASSLCSFFB.BoundaryChoices{DataStruct.ASSLCSFFB.BoundaryChoice}, 'Specify limits by time from beginning'))
                FFSyllBoundaries(Indices(i), 1) = ceil(DataStruct.ASSLCSFFB.StartLimit * Fs/1000);
                FFSyllBoundaries(Indices(i), 2) = length(Syllable) - ceil(DataStruct.ASSLCSFFB.EndLimit * Fs/1000);
            else
                if (strfind(DataStruct.ASSLCSFFB.BoundaryChoices{DataStruct.ASSLCSFFB.BoundaryChoice}, 'Specify limits by time from end'))
                    FFSyllBoundaries(Indices(i), 1) = length(Syllable) - ceil(DataStruct.ASSLCSFFB.StartLimit * Fs/1000);
                    FFSyllBoundaries(Indices(i), 2) = length(Syllable) - ceil(DataStruct.ASSLCSFFB.EndLimit * Fs/1000);
                end
            end
        end
    end
 
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

    if (i <= DataStruct.ASSLCSFFB.NumExamples)
        PlotSpectrogramInAxis_SongVar(Syllable, Time, Fs, DataStruct.axes1)
        axes(DataStruct.axes1);
        hold on;
        axis([0 1.1*Time(end) 300 8000]);
        
        plot(TempFF_x, TempFF, 'b', 'LineWidth', 2); 
        LeftBoundaryLine{i} = plot([(Time(1) + FFSyllBoundaries(Indices(i),1)/Fs) (Time(1) + FFSyllBoundaries(Indices(i),1)/Fs)], [300 8000], 'r--', 'LineWidth', 2);
        RightBoundaryLine{i} = plot([(Time(1) + FFSyllBoundaries(Indices(i),2)/Fs) (Time(1) + FFSyllBoundaries(Indices(i),2)/Fs)], [300 8000], 'g--', 'LineWidth', 2);
        plot(TempFF_x(StartIndex:EndIndex), FF(Indices(i)), 'g--', 'LineWidth', 2);
    end    
end

disp('Finished calculating boundaries');
