function [DataStruct] = ASSLReCalculateSyllFF(DataStruct, MaxFF, MinFF)

for i = 1:length(DataStruct.DataStruct.FileName),
    if (mod(i,10) == 0)
        fprintf('>%d', i);
    end
    [RawData, Fs] = ASSLGetRawData(DataStruct.DataStruct.DirName, DataStruct.DataStruct.FileName{i}, DataStruct.DataStruct.FileType, DataStruct.DataStruct.SongChanNo);
    % Parameters P for calculating FF for yin
    % sr for sample rate
    % maxf0 for max ff
    % minf0 for min ff
    P.sr = Fs;
    if (~isempty(MaxFF))
        P.maxf0 = MaxFF;
    end
    if (~isempty(MinFF))
        P.minf0 = MinFF;
    else
        P.minf0 = 300; % Default of 300 Hz for min FF
    end
    
    Temp_FF = yin(RawData, P);
    
    FF_T = linspace(1/Fs, length(RawData)/Fs, length(Temp_FF.f0));
    Temp_FF = 2.^Temp_FF.f0 * 440;
    
    Sylls = find(DataStruct.DataStruct.SyllLabels{i} == DataStruct.ASSLCSFFB.UniqueSyllLabels(DataStruct.ASSLCSFFB.SyllIndex));
    for j = 1:length(Sylls),
        StartIndex = find(FF_T <= DataStruct.DataStruct.SyllOnsets{i}(Sylls(j))/1000, 1, 'last');
        if (isempty(StartIndex))
            StartIndex = 1;
        end
    
        EndIndex = find(FF_T >= DataStruct.DataStruct.SyllOffsets{i}(Sylls(j))/1000, 1, 'first');
        if (isempty(EndIndex))
            EndIndex = length(FF_T);
        end
    
        DataStruct.DataStruct.Raw.FundamentalFrequency{i}{Sylls(j)} = Temp_FF(StartIndex:EndIndex);
    end
end        
fprintf('\n');