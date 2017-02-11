function [] = IntroNoteAmplitudePlotsAllMotif(Dir, IntroNoteResults, FileType, AlignmentPoint, PreTime, PostTime, ColorScale, BoutType, RepresentationType, ScaleForIntroNotes, varargin)

% Alignment Point is a number specifying the point where different bouts
% have to be aligned.
% 0 - refers to alignment on the onset of the first syllable of the motif
% -1 - refers to alignment on the onset of the last intro note
% -2 - refers to alignment on the onset of the second last intro
% -1000 - refers to alignment on the first intro note - in this case no
% trials will be dropped
% note
% and so on .....
%
% If the trial does not have the required number of intro notes for
% alignment, it will be excluded. For eg. if the alignment point is set at
% 'secondlast' and a bout has only one intro note then it will be excluded

% ======================================================================= %

if (nargin > 10)
    IntroNoteSylls = varargin{1};
end

switch (BoutType)
    
    case 'All'
        
        PresentDir = pwd;
        
        TotalBouts = 0;
        NewFs = 1000;
        AmpTime = -PreTime:1/NewFs:PostTime;
        AmpTime(end) = [];
        ZeroIndex = find(AmpTime == 0);
        if (isempty(ZeroIndex))
            ZeroIndex = length(AmpTime);
        end
        
        cd(Dir);
        ProgressBar = waitbar(0, 'Load up song files and generating colour representations');
        
        for i = 1:length(IntroNoteResults.NoofINs),
            TotalBouts = TotalBouts + 1;
            %disp(IntroNoteResults.BoutDetails(i).SongFile);
            
            if (strfind(FileType,'obs'))
                [Song,Fs] = soundin_copy(Dir, IntroNoteResults.BoutDetails(i).SongFile, 'obs0r');
                Song = Song * 5/32768;
            else
                if (strfind(FileType,'wav'));
                    [Song, Fs] = wavread(IntroNoteResults.BoutDetails(i).SongFile);
                else 
                    if (strfind(FileType, 'okrank'))
                        [Song, Fs] = ReadOKrankData(Dir, IntroNoteResults.BoutDetails(i).SongFile, 1);
                    end
                end
            end
            Time = [1:1:length(Song)]/Fs;
            
            AmplitudeFs = NewFs;
            AmplitudeTime = Time(1):1/AmplitudeFs:Time(end);
            if (size(AmplitudeTime,1) < size(AmplitudeTime,2))
                AmplitudeTime = AmplitudeTime';
            end
            Amplitude = zeros(size(AmplitudeTime));

            for j = 1:length(IntroNoteResults.BoutDetails(i).labels),
                if (ScaleForIntroNotes == 1)
                    if (IntroNoteResults.BoutDetails(i).labels(j) ~= 'i')
                        continue;
                    end
                end
                INStartIndex = find(AmplitudeTime <= (IntroNoteResults.BoutDetails(i).onsets(j)), 1, 'last');
                INEndIndex = find(AmplitudeTime <= (IntroNoteResults.BoutDetails(i).offsets(j)), 1, 'last');
                if (~isempty(strfind(RepresentationType, 'Syllable')))
                    INFlag = 0;
                    for IN = 1:length(IntroNoteSylls),
                        if (IntroNoteResults.BoutDetails(i).labels(j) == IntroNoteSylls(IN))
                            INFlag = 1;
                            break;
                        end
                    end
                    if (INFlag == 1)
                        Amplitude(INStartIndex:INEndIndex) = 1;
                    else
                        Amplitude(INStartIndex:INEndIndex) = 0.5;
                    end
                else
                    Amplitude(INStartIndex:INEndIndex) = IntroNoteResults.BoutDetails(i).Feats(j,RepresentationType);
                end
            end
            
            NoofINs(TotalBouts) = IntroNoteResults.NoofINs(i);
            if (NoofINs(TotalBouts) > 1)
                LastINInterval(TotalBouts) = IntroNoteResults.BoutDetails(i).onsets(IntroNoteResults.INs{i}(end)) - IntroNoteResults.BoutDetails(i).offsets(IntroNoteResults.INs{i}(end - 1));
            else
                LastINInterval(TotalBouts) = NaN;
            end
            
            AlignmentTime = IntroNoteResults.BoutDetails(i).onsets(IntroNoteResults.MotifStartIndex(i));
            AlignmentIndex = find(AmplitudeTime <= AlignmentTime, 1, 'last');
            StartIndex = AlignmentIndex - (ZeroIndex - 1);
            EndIndex = AlignmentIndex + (length(AmpTime) - ZeroIndex);
            
            BoutAmplitudes{TotalBouts} = [(AmplitudeTime(StartIndex:EndIndex) - AlignmentTime) Amplitude(StartIndex:EndIndex)];
            
            WithinBoutIndices = find(IntroNoteResults.WithinBoutINBoutIndices == i);
            for j = 1:length(WithinBoutIndices),
                TotalBouts = TotalBouts + 1;
                AlignmentTime = IntroNoteResults.BoutDetails(i).onsets(IntroNoteResults.WithinBoutNoofINs(WithinBoutIndices(j), 4));
                AlignmentIndex = find(AmplitudeTime <= AlignmentTime, 1, 'last');
                StartIndex = AlignmentIndex - (ZeroIndex - 1);
                EndIndex = AlignmentIndex + (length(AmpTime) - ZeroIndex);
             
                NoofINs(TotalBouts) = IntroNoteResults.WithinBoutNoofINs(WithinBoutIndices(j),1);
                if (NoofINs(TotalBouts) > 1)
                    LastINInterval(TotalBouts) = IntroNoteResults.BoutDetails(i).onsets(IntroNoteResults.WithinBoutINs{WithinBoutIndices(j)}(end)) - IntroNoteResults.BoutDetails(i).offsets(IntroNoteResults.WithinBoutINs{WithinBoutIndices(j)}(end - 1));
                else
                    LastINInterval(TotalBouts) = NaN;
                end
                BoutAmplitudes{TotalBouts} = [(AmplitudeTime(StartIndex:EndIndex) - AlignmentTime) Amplitude(StartIndex:EndIndex)];
            end
            waitbar(i/length(IntroNoteResults.NoofINs), ProgressBar, ['All bouts : ', IntroNoteResults.BoutDetails(i).SongFile]);
        end
        
        
        BoutNo = 0;
        for i = min(NoofINs):1:max(NoofINs),
            TrialIndices = find(NoofINs == i);
            TrialAmplitudes = BoutAmplitudes(TrialIndices);
            SortedIndices = 1:1:length(TrialIndices);
            for j = 1:length(SortedIndices),
                BoutNo = BoutNo + 1;
                AllMotifAmplitudes(BoutNo,:) = TrialAmplitudes{SortedIndices(j)}(:,2);
            end
        end
        
        figure
        if (ScaleForIntroNotes == 1)
            imagesc(AmpTime, [1:1:size(AllMotifAmplitudes,1)], (AllMotifAmplitudes), [min(AllMotifAmplitudes(find(AllMotifAmplitudes > 0))) max(AllMotifAmplitudes(find(AllMotifAmplitudes > 0)))]);
        else
            imagesc(AmpTime, [1:1:size(AllMotifAmplitudes,1)], (AllMotifAmplitudes));
        end
        set(gcf, 'Color', 'w');
        colormap(ColorScale);
        set(gca, 'XColor', 'k', 'YColor', 'k');
        xlabel('Time from the beginning of motif (sec)', 'FontSize', 12, 'FontName', 'Arial', 'Color', 'k');
        ylabel('Bout #', 'FontSize', 12, 'FontName', 'Arial', 'Color', 'k');
        set(gca, 'FontSize', 10, 'FontName', 'Arial', 'Color', 'k');
        set(gca, 'Box', 'off');
        disp('Finished Analysis');
        
        close(ProgressBar);
    
    case 'Beginning'
        
        PresentDir = pwd;

        TotalBouts = 0;
        NewFs = 1000;
        AmpTime = -PreTime:1/NewFs:PostTime;
        AmpTime(end) = [];
        ZeroIndex = find(AmpTime == 0);
        if (isempty(ZeroIndex))
            ZeroIndex = length(AmpTime);
        end
        
        clear BoutAmplitudes NoofINs LastINInterval;
        cd(Dir);

        ProgressBar = waitbar(0, 'Generating colour representation for bout beginnings');
        for i = 1:length(IntroNoteResults.NoofINs),
            TotalBouts = TotalBouts + 1;
            %disp(IntroNoteResults.BoutDetails(i).SongFile);

            if (strfind(FileType,'obs'))
                [Song,Fs] = soundin_copy(Dir, IntroNoteResults.BoutDetails(i).SongFile, 'obs0r');
                Song = Song * 5/32768;
            else
                if (strfind(FileType,'wav'));
                    [Song, Fs] = wavread(IntroNoteResults.BoutDetails(i).SongFile);
                else 
                    if (strfind(FileType, 'okrank'))
                        [Song, Fs] = ReadOKrankData(Dir, IntroNoteResults.BoutDetails(i).SongFile, 1);
                    end
                end
            end
            Time = [1:1:length(Song)]/Fs;

            AmplitudeFs = NewFs;
            AmplitudeTime = Time(1):1/AmplitudeFs:Time(end);
            if (size(AmplitudeTime,1) < size(AmplitudeTime,2))
                AmplitudeTime = AmplitudeTime';
            end
            Amplitude = zeros(size(AmplitudeTime));

            for j = 1:length(IntroNoteResults.BoutDetails(i).labels),
                INStartIndex = find(AmplitudeTime <= (IntroNoteResults.BoutDetails(i).onsets(j)), 1, 'last');
                INEndIndex = find(AmplitudeTime <= (IntroNoteResults.BoutDetails(i).offsets(j)), 1, 'last');
                if (~isempty(strfind(RepresentationType, 'Syllable')))
                    INFlag = 0;
                    for IN = 1:length(IntroNoteSylls),
                        if (IntroNoteResults.BoutDetails(i).labels(j) == IntroNoteSylls(IN))
                            INFlag = 1;
                            break;
                        end
                    end
                    if (INFlag == 1)
                        Amplitude(INStartIndex:INEndIndex) = 1;
                    else
                        Amplitude(INStartIndex:INEndIndex) = 0.5;
                    end
                else
                    Amplitude(INStartIndex:INEndIndex) = IntroNoteResults.BoutDetails(i).Feats(j,RepresentationType);
                end
            end
            
            NoofINs(TotalBouts) = IntroNoteResults.NoofINs(i);
            if (NoofINs(TotalBouts) > 1)
                LastINInterval(TotalBouts) = IntroNoteResults.BoutDetails(i).onsets(IntroNoteResults.INs{i}(end)) - IntroNoteResults.BoutDetails(i).offsets(IntroNoteResults.INs{i}(end - 1));
            else
                LastINInterval(TotalBouts) = NaN;
            end

            AlignmentTime = IntroNoteResults.BoutDetails(i).onsets(IntroNoteResults.MotifStartIndex(i));
            AlignmentIndex = find(AmplitudeTime <= AlignmentTime, 1, 'last');
            StartIndex = AlignmentIndex - (ZeroIndex - 1);
            if (StartIndex <= 0)
                PointsToBeAdded = abs(StartIndex) + 1;
                StartIndex = 1;
            else
                PointsToBeAdded = 0;
            end
            EndIndex = AlignmentIndex + (length(AmpTime) - ZeroIndex);

            if (PointsToBeAdded == 0)
                BoutAmplitudes{TotalBouts} = [(AmplitudeTime(StartIndex:EndIndex) - AlignmentTime) Amplitude(StartIndex:EndIndex)];
            else
                AmplitudeTimeFs = 1/(AmplitudeTime(StartIndex + 1) - AmplitudeTime(StartIndex));
                BoutAmplitudes{TotalBouts} = [[((AmplitudeTime(StartIndex) - (1:1:PointsToBeAdded)'/AmplitudeTimeFs) - AlignmentTime); (AmplitudeTime(StartIndex:EndIndex) - AlignmentTime)] [zeros(PointsToBeAdded, 1); Amplitude(StartIndex:EndIndex)]];
            end
            waitbar(i/length(IntroNoteResults.NoofINs), ProgressBar, ['Bout beginning : ', IntroNoteResults.BoutDetails(i).SongFile]);
        end

        clear TrialAmplitudes AllMotifAmplitudes;
        BoutNo = 0;
        for i = min(NoofINs):1:max(NoofINs),
            TrialIndices = find(NoofINs == i);
            TrialAmplitudes = BoutAmplitudes(TrialIndices);
            %if (i > 1)
            %    TrialLastInterval = LastINInterval(TrialIndices);
            %    [SortedIntervals, SortedIndices] = sort(TrialLastInterval);
            %else
            SortedIndices = 1:1:length(TrialIndices);
            %end
            for j = 1:length(SortedIndices),
                BoutNo = BoutNo + 1;
                AllMotifAmplitudes(BoutNo,:) = TrialAmplitudes{SortedIndices(j)}(:,2);
            end
        end

        figure
        if (ScaleForIntroNotes == 1)
            imagesc(AmpTime, [1:1:size(AllMotifAmplitudes,1)], (AllMotifAmplitudes), [min(AllMotifAmplitudes(find(AllMotifAmplitudes > 0))) max(AllMotifAmplitudes(find(AllMotifAmplitudes > 0)))]);
        else
            imagesc(AmpTime, [1:1:size(AllMotifAmplitudes,1)], (AllMotifAmplitudes));
        end
        set(gcf, 'Color', 'w');
        colormap(ColorScale);
        set(gca, 'XColor', 'k', 'YColor', 'k');
        xlabel('Time from the beginning of motif (sec)', 'FontSize', 12, 'FontName', 'Arial', 'Color', 'k');
        ylabel('Bout #', 'FontSize', 12, 'FontName', 'Arial', 'Color', 'k');
        set(gca, 'FontSize', 10, 'FontName', 'Arial', 'Color', 'k');
        set(gca, 'Box', 'off');
        set(gca, 'YColor', 'w');
        disp('Finished Analysis');
        close(ProgressBar);

    case 'Within'
        
        PresentDir = pwd;

        TotalBouts = 0;
        NewFs = 1000;
        AmpTime = -PreTime:1/NewFs:PostTime;
        AmpTime(end) = [];
        ZeroIndex = find(AmpTime == 0);
        
        if (isempty(ZeroIndex))
            ZeroIndex = length(AmpTime);
        end
        
        clear BoutAmplitudes NoofINs LastINInterval;

        cd(Dir);
        ProgressBar = waitbar(0, 'Generating colour representations for within bouts');
        for i = 1:length(IntroNoteResults.NoofINs),
            %TotalBouts = TotalBouts + 1;
            %disp(IntroNoteResults.BoutDetails(i).SongFile);

            if (strfind(FileType,'obs'))
                [Song,Fs] = soundin_copy(Dir, IntroNoteResults.BoutDetails(i).SongFile, 'obs0r');
                Song = Song * 5/32768;
            else
                if (strfind(FileType,'wav'));
                    [Song, Fs] = wavread(IntroNoteResults.BoutDetails(i).SongFile);
                else 
                    if (strfind(FileType, 'okrank'))
                        [Song, Fs] = ReadOKrankData(Dir, IntroNoteResults.BoutDetails(i).SongFile, 1);
                    end
                end
            end
            Time = [1:1:length(Song)]/Fs;

            AmplitudeFs = NewFs;
            AmplitudeTime = Time(1):1/AmplitudeFs:Time(end);
            if (size(AmplitudeTime,1) < size(AmplitudeTime,2))
                AmplitudeTime = AmplitudeTime';
            end
            Amplitude = zeros(size(AmplitudeTime));

            for j = 1:length(IntroNoteResults.BoutDetails(i).labels),
                INStartIndex = find(AmplitudeTime <= (IntroNoteResults.BoutDetails(i).onsets(j)), 1, 'last');
                INEndIndex = find(AmplitudeTime <= (IntroNoteResults.BoutDetails(i).offsets(j)), 1, 'last');
                if (~isempty(strfind(RepresentationType, 'Syllable')))
                    INFlag = 0;
                    for IN = 1:length(IntroNoteSylls),
                        if (IntroNoteResults.BoutDetails(i).labels(j) == IntroNoteSylls(IN))
                            INFlag = 1;
                            break;
                        end
                    end
                    if (INFlag == 1)
                        Amplitude(INStartIndex:INEndIndex) = 1;
                    else
                        Amplitude(INStartIndex:INEndIndex) = 0.5;
                    end
                else
                    Amplitude(INStartIndex:INEndIndex) = IntroNoteResults.BoutDetails(i).Feats(j,RepresentationType);
                end
            end

            WithinBoutIndices = find(IntroNoteResults.WithinBoutINBoutIndices == i);
            for j = 1:length(WithinBoutIndices),
                TotalBouts = TotalBouts + 1;
                if (AlignmentPoint == 0)
                    AlignmentTime = IntroNoteResults.BoutDetails(i).onsets(IntroNoteResults.WithinBoutNoofINs(WithinBoutIndices(j), 4));
                else
                    if (AlignmentPoint == -1000)
                        if (isempty(IntroNoteResults.WithinBoutINs{WithinBoutIndices(j)}))
                            AlignmentTime = IntroNoteResults.BoutDetails(i).onsets(IntroNoteResults.WithinBoutNoofINs(WithinBoutIndices(j), 4));
                        else
                            AlignmentTime = IntroNoteResults.BoutDetails(i).onsets(IntroNoteResults.WithinBoutINs{WithinBoutIndices(j)}(1));
                        end
                    else
                        if (AlignmentPoint == -2000)
                            AlignmentTime = IntroNoteResults.BoutDetails(i).offsets(IntroNoteResults.WithinBoutNoofINs(WithinBoutIndices(j), 7));
                        end
                    end
                end
                AlignmentIndex = find(AmplitudeTime <= AlignmentTime, 1, 'last');
                StartIndex = AlignmentIndex - (ZeroIndex - 1);
                EndIndex = AlignmentIndex + (length(AmpTime) - ZeroIndex);

                NoofINs(TotalBouts) = IntroNoteResults.WithinBoutNoofINs(WithinBoutIndices(j),1);
                %if (NoofINs(TotalBouts) > 1)
                %    LastINInterval(TotalBouts) = IntroNoteResults.BoutDetails(i).onsets(IntroNoteResults.WithinBoutINs{WithinBoutIndices(j)}(end)) - IntroNoteResults.BoutDetails(i).offsets(IntroNoteResults.WithinBoutINs{WithinBoutIndices(j)}(end - 1));
                %else
                %    LastINInterval(TotalBouts) = NaN;
                %end
                LastINInterval(TotalBouts) = IntroNoteResults.WithinBoutNoofINs(WithinBoutIndices(j), 5);
                BoutAmplitudes{TotalBouts} = [(AmplitudeTime(StartIndex:EndIndex) - AlignmentTime) Amplitude(StartIndex:EndIndex)];
            end
            waitbar(i/length(IntroNoteResults.NoofINs), ProgressBar, ['Within bouts : ', IntroNoteResults.BoutDetails(i).SongFile]);
        end

        clear TrialAmplitudes AllMotifAmplitudes;
        BoutNo = 0;
        for i = min(NoofINs):1:max(NoofINs),
            TrialIndices = find(NoofINs == i);
            TrialAmplitudes = BoutAmplitudes(TrialIndices);
            %if (i > 1)
                TrialLastInterval = LastINInterval(TrialIndices);
                [SortedIntervals, SortedIndices] = sort(TrialLastInterval);
            %else
            %    SortedIndices = 1:1:length(TrialIndices);
            %end
            for j = 1:length(SortedIndices),
                BoutNo = BoutNo + 1;
                AllMotifAmplitudes(BoutNo,:) = TrialAmplitudes{SortedIndices(j)}(:,2);
            end
        end

        figure
        imagesc(AmpTime, [1:1:size(AllMotifAmplitudes,1)], (AllMotifAmplitudes));
        set(gcf, 'Color', 'w');
        colormap(ColorScale);
        set(gca, 'XColor', 'k', 'YColor', 'k');
        xlabel('Time from the beginning of motif (sec)', 'FontSize', 12, 'FontName', 'Arial', 'Color', 'k');
        ylabel('Bout #', 'FontSize', 12, 'FontName', 'Arial', 'Color', 'k');
        set(gca, 'FontSize', 10, 'FontName', 'Arial', 'Color', 'k');
        set(gca, 'Box', 'off');
        set(gca, 'YColor', 'w');
        close(ProgressBar);
end

disp('Finished Analysis');
