function [] = ParseBoutStatisticsFiles(BirdName, Date)

PresentDir = pwd;

if (ispc)
    cd('E:\RaghavData\BoutStatistics');
else
    cd('/home/raghav/BoutStatistics');
end

BoutStatsFiles = dir('*.log');
DMatchIndex = 1;
UMatchIndex = 1;
for i = 1:length(BoutStatsFiles),
    Temp = textread([BoutStatsFiles(i).name], '%s', 'delimiter', '\n');
    if (~isempty(strfind(Temp{1}, BirdName)) && ~isempty(strfind(Temp{2}, Date)))
        clear MaxMatches SongFileNames;
        MaxMatchLines = find(cellfun(@length, strfind(Temp, 'Maximum')));
        MatchValIndex = 1;
        for j = 1:length(MaxMatchLines),
            [s, s1, s2, s3, Max1, s4, s5, Max2] = strread(Temp{MaxMatchLines(j)}, '%s %s %s %s %f %s %s %f');
            MaxMatches(MatchValIndex) = Max1;
            MatchValIndex = MatchValIndex + 1;
            clear Max1;
        end
        if (exist('MaxMatches', 'var'))
            if (isempty(MaxMatches))
                MaxMatches = -100;
            end
        else
            MaxMatches = -100;
        end
        if (~isempty(strfind(Temp{2}, 'alone')) || ~isempty(strfind(Temp{2}, 'undir')) || ~isempty(strfind(Temp{2}, 'Undir')) || ~isempty(strfind(Temp{2}, 'UnDir')))
            UMatches{UMatchIndex} = BoutStatsFiles(i).name;
            UMatchIndex = UMatchIndex + 1;
            disp(['Undir song: ', BoutStatsFiles(i).name, ' - Median match is ', num2str(median(MaxMatches))]);
        else
            DMatches{DMatchIndex} = BoutStatsFiles(i).name;
            DMatchIndex = DMatchIndex + 1;
            disp(['Dir song: ', BoutStatsFiles(i).name, ' - Median match is ', num2str(median(MaxMatches))]);
        end
    end
    clear Temp;
end

cd(PresentDir);