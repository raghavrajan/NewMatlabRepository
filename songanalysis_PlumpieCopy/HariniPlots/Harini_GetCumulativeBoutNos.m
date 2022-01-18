function [] = Harini_GetCumulativeBoutNos(CSVTextFile)

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
    end
    [SortedVals, SortedIndices] = sort(DateNumber{i});
    Date{i} = Date{i}(SortedIndices);
    DateNumber{i} = DateNumber{i}(SortedIndices);
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

for i = 1:length(UniqueBirds),
    figure;
    set(gcf, 'Color', 'w');
    set(gcf, 'Position', [318 66 1700 750]);
    p = panel();
    p.pack({2/3 1/3});
    p.fontsize = 10;
    p.fontweight = 'bold';
    p(1).select();
    Temp = [];
    for j = 1:length(DateNumber{i}),
        Temp(:,j) = (DateNumber{i}(j) - DateNumber{i});
    end
    imagesc(Temp);
    set(gca, 'XTick', 1:1:length(DateNumber{i}), 'XTickLabel', Date{i}, 'XTickLabelRotation', 30);
    set(gca, 'YTick', 1:1:length(DateNumber{i}), 'YTickLabel', Date{i}, 'YTickLabelRotation', 0);
    %caxis([-26 26]);
    %colormap([[0 0 1]; repmat([0.8 0.8 0.8], 51, 1); [0 0 1]]);
    colorbar
    xlabel('Recording days');
    ylabel('Recording days');
    title([UniqueBirds{i}, ': time difference between successive recording sessions (days)']);
    axis tight;
    
    p(2).select();
    plot(cumsum(NumSongBouts{i}), 'ko-');
    axis tight;
    Temp = axis;
    Temp = [0.5 (length(NumSongBouts{i})+0.5) 0 1.02*Temp(4)];
    axis(Temp);
    set(gca, 'XTick', 1:1:length(NumSongBouts{i}), 'XTickLabel', Date{i}, 'XTickLabelRotation', 30);
    xlabel('Recording days');
    ylabel('Cumulative # of song bouts');
    title([UniqueBirds{i}, ': cumulative # of song bouts across successive recording sessions']);
    
    p.de.margin = 25;
    p.marginleft = 25;
    p.margintop = 10;
    p.marginbottom = 20;
    OutputDir = '/home/raghav/StudentRelated/Harini';
    set(gcf, 'PaperPositionMode', 'auto');
    print(fullfile(OutputDir, [UniqueBirds{i}, '.RecordingSessionTimeDifferences.CumulativeNumSongBouts.png']), '-dpng');
end

disp('Got bout nos');


