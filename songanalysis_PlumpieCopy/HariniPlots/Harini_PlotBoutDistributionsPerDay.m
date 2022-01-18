function [] = Harini_PlotBoutDistributionsPerDay(CSVTextFile)

ContinuousFileTime = 30; % in seconds
InterBoutInterval = 2; % in seconds

% First pre-processing step
[BirdParameters, Flag] = Harini_ProcessSongData(CSVTextFile, 2000);

% Now I need to get the first valid song bout and calculate the latency to
% the first song

% I also need the total number of song bouts. In this case, there are song
% bouts that sometimes spread over two files. In this case I need to check
% the time of the file and see whether two files are continuous or not.

% For this, I need to first get the times for each of the files
for i = 1:length(BirdParameters),
    for j = 1:length(BirdParameters(i).SongFileNames),
        ExtensionIndex = strfind(BirdParameters(i).SongFileNames{j}, '.wav');
        FileTimeString = BirdParameters(i).SongFileNames{j}((ExtensionIndex - 6):(ExtensionIndex - 1)); % string from filename in the format hhmmss (h for hour, m for minutes, s for seconds)
        BirdParameters(i).FileTime(j) = (str2double(FileTimeString(1:2))) + (str2double(FileTimeString(3:4))/60) + (str2double(FileTimeString(5:6))/3600); % in hours
    end
end

% Get unique BirdNames to get an idea of how many birds there are
for i = 1:length(BirdParameters),
    BirdNames{i} = BirdParameters(i).BirdName;
    MicrophoneTypes{i} = BirdParameters(i).Microphone;
end

UniqueBirds = unique(BirdNames);

% Now to get the number of song bouts in total for each bird for each day
for i = 1:length(UniqueBirds),
    Matches = strmatch(UniqueBirds{i}, BirdNames);
    RecordingDates = [];
    for j = 1:length(Matches),
        Date{i}{j} = BirdParameters(Matches(j)).DataLabel;
        DateNumber{i}(j) = datenum(Date{i}{j}, 'ddmmyy');
        Condition{i}{j} = BirdParameters(Matches(j)).Condition;
        Microphone{i}{j} = BirdParameters(Matches(j)).Microphone;
    end
    [SortedVals, SortedIndices] = sort(DateNumber{i});
    Date{i} = Date{i}(SortedIndices);
    DateNumber{i} = DateNumber{i}(SortedIndices);
    Condition{i} = Condition{i}(SortedIndices);
    Microphone{i} = Microphone{i}(SortedIndices);
    
    for j = SortedIndices(:)',
        SkipNextBout = 0;
        ValidSongBouts = find(BirdParameters(Matches(j)).Bouts(:,7) == 1);
        NumSongBouts{i}(j) = 0;
        BirdParameters(Matches(j)).ValidSongBouts = [];
        if (isempty(ValidSongBouts))
            continue;
        end
        
        if (length(ValidSongBouts) == 1)
            BirdParameters(Matches(j)).ValidSongBouts = ValidSongBouts;
            NumSongBouts{i}(j) = 1;
        else
            for k = ValidSongBouts(:)',
                if (SkipNextBout == 1)
                    SkipNextBout = 0;
                    continue;
                end
                
                SongFileIndex = BirdParameters(Matches(j)).Bouts(k,3);
                BoutOnsetTime = BirdParameters(Matches(j)).Bouts(k,5)/1000; % in seconds
                BoutOffsetTime = BirdParameters(Matches(j)).Bouts(k,6)/1000; % in seconds
                
                if ((BoutOnsetTime >= InterBoutInterval) && ((ContinuousFileTime - BoutOffsetTime) >= InterBoutInterval))
                    NumSongBouts{i}(j) = NumSongBouts{i}(j) + 1;
                    BirdParameters(Matches(j)).ValidSongBouts(end+1) = k;
                else
                    if ((ContinuousFileTime - BoutOffsetTime) >= InterBoutInterval)
                        NumSongBouts{i}(j) = NumSongBouts{i}(j) + 1;
                        BirdParameters(Matches(j)).ValidSongBouts(end+1) = k;
                    else
                        % First check whether there is a next file. If not,
                        % just include the bout and continue
                        if (k == ValidSongBouts(end))
                            NumSongBouts{i}(j) = NumSongBouts{i}(j) + 1;
                            BirdParameters(Matches(j)).ValidSongBouts(end+1) = k;
                            continue;
                        end
                        
                        % Now to check whether next file in the
                        % chronological order is the next file based on
                        % time too. If no, then count as a bout.
                        % Now to check if next valid song bout is in
                        % the next file and is close to the start
                        CurrentBout = find(ValidSongBouts == k);
                        NextBout = ValidSongBouts(CurrentBout + 1);
                        if ((BirdParameters(Matches(j)).Bouts(NextBout,3) - BirdParameters(Matches(j)).Bouts(k,3)) == 1)
                            % Now if the time before the start of the
                            % next bout and the end of the current bout
                            % is less than 2s, then merge the two bouts
                            CurrentBoutOffsetTime = BirdParameters(Matches(j)).FileTime(SongFileIndex) + BoutOffsetTime/3600;
                            NextBoutOnsetTime = BirdParameters(Matches(j)).FileTime(BirdParameters(Matches(j)).Bouts(NextBout,3)) + BirdParameters(Matches(j)).Bouts(NextBout,5)/(1000*3600);
                            if ((NextBoutOnsetTime - CurrentBoutOffsetTime)*3600 < InterBoutInterval)
                                NumSongBouts{i}(j) = NumSongBouts{i}(j) + 1;
                                BirdParameters(Matches(j)).ValidSongBouts(end+1) = k;
                                SkipNextBout = 1;
                            else
                                NumSongBouts{i}(j) = NumSongBouts{i}(j) + 1;
                                BirdParameters(Matches(j)).ValidSongBouts(end+1) = k;
                            end
                        else
                            NumSongBouts{i}(j) = NumSongBouts{i}(j) + 1;
                            BirdParameters(Matches(j)).ValidSongBouts(end+1) = k;
                        end
                    end
                end
            end
        end
    end
end

Colours = 'rgbcmk';

AllConditions = {'L0' 'L1' 'L2' 'L3' 'L4' 'UN'};

for i = 1:length(UniqueBirds),
    FirstConditions = zeros(6,1); % Just a counter to keep track of the first of a condition and use that for the legend
    
    figure;
    set(gcf, 'Color', 'w');
    set(gcf, 'Position', [95 414 1700 500]);
    p = panel();
    p.pack({1});
    p.fontsize = 12;
    p.fontweight = 'bold';
    p(1).select();
    hold on;
    
    clear DateLabelString
    [UniqueDates, UniqueDateIndices] = unique(DateNumber{i});
    UniqueDates = Date{i}(UniqueDateIndices);
    for j = 1:length(UniqueDates),
        Matches = strmatch(UniqueDates{j}, Date{i});
        if (j == 1)
            FirstDayNumber = DateNumber{i}(Matches(1));
            DateLabelString{j} = ['Day #1: ', Date{i}{Matches(1)}];
        else
            DateLabelString{j} = ['Day #', num2str(1 + (DateNumber{i}(Matches(1)) - FirstDayNumber)), ': ', Date{i}{Matches(1)}];
        end
        DateLabelString{j} = [DateLabelString{j}, ': ', Microphone{i}{Matches(1)}];
                
        for k = 1:length(Matches),
            switch (Condition{i}{Matches(k)})
                case 'L0'
                    ColourIndex = 1;
                case 'L1'
                    ColourIndex = 2;
                case 'L2'
                    ColourIndex = 3;
                case 'L3'
                    ColourIndex = 4;
                case 'L4'
                    ColourIndex = 5;
                case 'UN'
                    ColourIndex = 6;
            end
            
            if (FirstConditions(ColourIndex) == 0)
                LegendBar{i}(ColourIndex) = bar(j+((k-1)/5), NumSongBouts{i}(Matches(k)), 0.15);
                set(LegendBar{i}(ColourIndex), 'FaceColor', Colours(ColourIndex), 'EdgeColor', Colours(ColourIndex));
                FirstConditions(ColourIndex) = 1;
            else
                TempBar = bar(j+((k-1)/5), NumSongBouts{i}(Matches(k)), 0.15);
                set(TempBar, 'FaceColor', Colours(ColourIndex), 'EdgeColor', Colours(ColourIndex));
            end
        end
    end
    
    ValidConditions = find(FirstConditions);
    legend(LegendBar{i}(ValidConditions), AllConditions(ValidConditions));
    axis tight;
    Temp = axis;
    Temp = [0.5 (length(UniqueDates)+0.5) 0 1.02*Temp(4)];
    axis(Temp);
    set(gca, 'XTick', 1:1:length(UniqueDates), 'XTickLabel', DateLabelString, 'XTickLabelRotation', 30);
    xlabel('Recording days');
    ylabel('# of song bouts');
    title([UniqueBirds{i}, ': # of song bouts across successive recording sessions']);
    
    p.de.margin = 25;
    p.marginleft = 25;
    p.margintop = 10;
    p.marginbottom = 50;
    OutputDir = '/home/raghav/StudentRelated/Harini';
    set(gcf, 'PaperPositionMode', 'auto');
    print(fullfile(OutputDir, [UniqueBirds{i}, '.RecordingSessionDateCondition.NumSongBouts.png']), '-dpng');
end

disp('Got bout nos');


