function [SyllFF] = ASSLCalculateFFWithOnsets(Song, Time, Fs, Onsets, Offsets, StartPercentage, EndPercentage, MaxFF, MinFF)

if (~isempty(Onsets))
    % Parameters P for calculating FF for yin
    % sr for sample rate
    % maxf0 for max ff
    % minf0 for min ff
    P.sr = Fs;
    P.minf0 = MinFF; 
    P.maxf0 = MaxFF;
    
    % A dirty fix for some error I was getting - am not sure why the error
    % was coming so used a try catch routine to bypass this error for the n
    try
        FF = yin(Song, P);
        FF_T = linspace(Time(1), Time(end), length(FF.f0)); % in sec
        FF = 2.^FF.f0 * 440;
    catch
        FF = ones(size(Onsets)) * NaN;
    end
end

if (isempty(Onsets))
    SyllFF = [];
end

for i = 1:length(Onsets),
    % For each of these syllables, I have use the start percentage and end
    % percentage as the times of starting and ending for calculating the
    % mean ff
    
    SyllLength = Offsets(i) - Onsets(i);
    StartTime = Onsets(i) + ((StartPercentage * SyllLength)/100);
    EndTime = Onsets(i) + ((EndPercentage * SyllLength)/100);
    
    StartIndex = find(FF_T <= StartTime, 1, 'last');
    if (isempty(StartIndex))
        StartIndex = 1;
    end
    
    EndIndex = find(FF_T >= EndTime, 1, 'first');
    if (isempty(EndIndex))
        EndIndex = length(FF_T);
    end
    
    if (length(find(isnan(FF))) ~= length(FF))
        SyllFF(i,1) = nanmedian(FF(StartIndex:EndIndex)); % taking median since FF calculation is not entirely reliable and FF values jump around and so median will ensure that the erroneous values will not affect our FF value for each syllable
        SyllFF(i,2) = nanmean(FF(StartIndex:EndIndex));
    else
        SyllFF(i,1) = NaN;
        SyllFF(i,2) = NaN;
    end
end        
