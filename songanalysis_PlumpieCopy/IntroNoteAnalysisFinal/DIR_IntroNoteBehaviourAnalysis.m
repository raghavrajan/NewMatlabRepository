function [IntroNoteResults] = DIR_IntroNoteBehaviourAnalysis(DataDir, RecFileDir, SongFileList, NoteFileDir, FileType, MotifSylls, InterBoutInterval, StartSylls, ContinuousOrDiscrete, ContinuousFileTime, varargin)

% Get bout details with a function ScreenBouts %

if (nargin > 10)
    ContinuousFiles = varargin{1};
    [BoutDetails, AllBoutDetails] = DIR_ScreenBouts(DataDir, RecFileDir, FileType, NoteFileDir, SongFileList, InterBoutInterval, MotifSylls, ContinuousOrDiscrete, ContinuousFileTime, ContinuousFiles);
else
    [BoutDetails, AllBoutDetails] = DIR_ScreenBouts(DataDir, RecFileDir, FileType, NoteFileDir, SongFileList, InterBoutInterval, MotifSylls, ContinuousOrDiscrete, ContinuousFileTime);
end

% Now, find the first motif syllable in a bout and see if it is preceded by
% intro notes or not - get the proportion of times the first motif syllable
% is preceded by IN

INBouts = 0;
IntroNoteResults.WithinBoutNoofINs = [];
INWithinBoutNo = 0;

for i = 1:length(BoutDetails),
    FirstSyll(i) = BoutDetails(i).labels(1);
    Matches = [];
    for j = 1:length(MotifSylls),
        TempMatches = find((BoutDetails(i).labels == MotifSylls(j)) | (BoutDetails(i).labels == upper(MotifSylls(j))));
        if (~isempty(TempMatches))
            Matches(j) = TempMatches(1);
        else
            Matches(j) = -100;
        end
    end
    Matches = min(Matches(find(Matches > 0)));
    if (Matches > 1)
        IntroNoteResults.NoofSylls(i) = Matches - 1;
        IN_Indices = find(BoutDetails(i).labels(1:Matches) == 'i');
        if (isempty(find(diff(IN_Indices) > 1)))
            IntroNoteResults.INs{i} = IN_Indices;
        else
            DiffIndex = find(diff(IN_Indices) > 1);
            IntroNoteResults.INs{i} = IN_Indices((DiffIndex(end) + 1):end);
        end
        IntroNoteResults.MotifStartIndex(i) = Matches;
        IntroNoteResults.NoofINs(i) = length(IntroNoteResults.INs{i});
        IntroNoteResults.INFeats{i} = BoutDetails(i).Feats(IntroNoteResults.INs{i},:);
        IntroNoteResults.Sylls{i} = 1:1:Matches-1;
        IntroNoteResults.SyllFeats{i} = BoutDetails(i).Feats(IntroNoteResults.Sylls{i},:);
        IntroNoteResults.TimeToMotif(i) = BoutDetails(i).onsets(Matches) - BoutDetails(i).onsets(1);
        SyllBeforeMotif(i) = BoutDetails(i).labels(Matches -1);
    else
        if (Matches == 1)
            IntroNoteResults.MotifStartIndex(i) = 1;
            IntroNoteResults.NoofINs(i) = 0;
            IntroNoteResults.NoofSylls(i) = 0;
            IntroNoteResults.TimeToMotif(i) = BoutDetails(i).onsets(Matches) - BoutDetails(i).onsets(1);
            SyllBeforeMotif(i) = 'Q';
        end
    end
    
    Matches = [];
    for j = 1:length(StartSylls),
        Matches = [Matches find(BoutDetails(i).labels == StartSylls(j))];
    end
    Matches = sort(Matches);
    Matches(find(diff(Matches) == 1) + 1) = [];
    IntroNoteResults.BoutNoofMotifs(i) = length(Matches);
    for j = 2:length(Matches),
        MotifIntroSyllsFlag = [];
        MotifSyllsFlag = [];
        for k = Matches(j-1):Matches(j),
            if (isempty(strfind([MotifSylls, 'i'], BoutDetails(i).labels(k))))
                MotifIntroSyllsFlag = [MotifIntroSyllsFlag 0];
            else
                MotifIntroSyllsFlag = [MotifIntroSyllsFlag 1];
            end
            if (isempty(strfind(MotifSylls, BoutDetails(i).labels(k))))
                MotifSyllsFlag = [MotifSyllsFlag 0];
            else
                MotifSyllsFlag = [MotifSyllsFlag 1];
            end
        end
        LastMotifSyll = Matches(j-1) + find(MotifSyllsFlag(1:end-1) == 1, 1, 'last') - 1;

        IntroNotes = find(BoutDetails(i).labels(Matches(j-1):Matches(j)) == 'i');
        % just verifying that the last intro note is just before the new
        % motif, otherwise these intro notes are not considered
        if (~isempty(IntroNotes))
            if ((IntroNotes(end) + Matches(j-1) - 1) ~= (Matches(j) - 1))
                IntroNotes = [];
            end
        end
        INWithinBoutNo = INWithinBoutNo + 1;
        IntroNoteResults.WithinBoutINBoutIndices(INWithinBoutNo) = i;

        Calls = find(MotifIntroSyllsFlag == 0);
        
        if (isempty(IntroNotes))
            [LongestInterval, LongestIntervalIndex] = max(BoutDetails(i).onsets([LastMotifSyll+1:Matches(j)]) - BoutDetails(i).offsets([LastMotifSyll:Matches(j)-1]));
            LongestIntervalIndex = LongestIntervalIndex + LastMotifSyll - 1;
            IntroNoteResults.WithinBoutNoofINs(INWithinBoutNo,:) = [0 min(MotifIntroSyllsFlag) (BoutDetails(i).onsets(Matches(j)) - BoutDetails(i).offsets(LastMotifSyll)) Matches(j) LongestInterval LastMotifSyll LongestIntervalIndex length(BoutDetails(i).labels(LongestIntervalIndex+1:Matches(j)-1))];
        else
            IntroNotes = IntroNotes + Matches(j-1) - 1;
            if (~isempty(find(diff(IntroNotes) > 1)))
                DiffIndex = find(diff(IntroNotes) > 1);
                IntroNotes = IntroNotes((DiffIndex(end) + 1):end);
            end
            [LongestInterval, LongestIntervalIndex] = max(BoutDetails(i).onsets([LastMotifSyll+1:IntroNotes(1)]) - BoutDetails(i).offsets([LastMotifSyll:IntroNotes(1)-1]));            
            LongestIntervalIndex = LongestIntervalIndex + LastMotifSyll - 1;
            IntroNoteResults.WithinBoutINs{INWithinBoutNo} = IntroNotes;
            IntroNoteResults.WithinBoutNoofINs(INWithinBoutNo,:) = [length(IntroNotes) min(MotifIntroSyllsFlag) (BoutDetails(i).onsets(IntroNotes(1)) - BoutDetails(i).offsets(LastMotifSyll)) Matches(j) LongestInterval LastMotifSyll LongestIntervalIndex length(BoutDetails(i).labels(LongestIntervalIndex+1:Matches(j)-1))];
            IntroNoteResults.WithinBoutINFeats{INWithinBoutNo} = BoutDetails(i).Feats(IntroNotes,:);
        end
        if (~isempty(length(BoutDetails(i).labels(LongestIntervalIndex+1:Matches(j)-1))))
            IntroNoteResults.WithinBoutSylls{INWithinBoutNo} = LongestIntervalIndex+1:Matches(j)-1;
            IntroNoteResults.WithinBoutSyllFeats{INWithinBoutNo} = BoutDetails(i).Feats(IntroNoteResults.WithinBoutSylls{INWithinBoutNo});
            IntroNoteResults.WithinBoutNoofSylls(INWithinBoutNo) = length(IntroNoteResults.WithinBoutSylls{INWithinBoutNo});
        end
    end
end

IntroNoteResults.FirstSyll = FirstSyll;
IntroNoteResults.SyllBeforeMotifs = SyllBeforeMotif;

disp('First Syllables of a bout and their proportions are as follows: ');
UniqueFirstSylls = unique(FirstSyll);
for i = 1:length(UniqueFirstSylls),
    disp([UniqueFirstSylls(i), ': ', num2str(length(find(FirstSyll == UniqueFirstSylls(i)))*100/length(FirstSyll)), '%']);
    IntroNoteResults.FirstSylls(i) = UniqueFirstSylls(i);
    IntroNoteResults.FirstSyllsProp(i) = length(find(FirstSyll == UniqueFirstSylls(i)))*100/length(FirstSyll);
end

disp('Syllables before the first motif syllable and their proportions are as follows: ');
UniqueSyllBeforeMotif = unique(SyllBeforeMotif);
for i = 1:length(UniqueSyllBeforeMotif),
    disp([UniqueSyllBeforeMotif(i), ': ', num2str(length(find(SyllBeforeMotif == UniqueSyllBeforeMotif(i)))*100/length(SyllBeforeMotif)), '%']);
    IntroNoteResults.SyllBeforeMotif(i) = UniqueSyllBeforeMotif(i);
    IntroNoteResults.SyllBeforeMotifProp(i) = length(find(SyllBeforeMotif == UniqueSyllBeforeMotif(i)))*100/length(SyllBeforeMotif);
end

% figure;
% Edges = 0:1:max(IntroNoteResults.NoofINs) + 1;
% plot(Edges, histc(NoofINs(find(FirstSyll == 'i')), Edges)/sum(histc(NoofINs(find(FirstSyll == 'i')), Edges)), 'rs-');
% hold on;
% plot(Edges, histc(NoofINs(find(FirstSyll ~= 'i')), Edges)/sum(histc(NoofINs(find(FirstSyll ~= 'i')), Edges)), 'ks-');

IntroNoteResults.MeanIN_INBouts = mean(IntroNoteResults.NoofINs(find(FirstSyll == 'i')));
IntroNoteResults.STDIN_INBouts = std(IntroNoteResults.NoofINs(find(FirstSyll == 'i')));
IntroNoteResults.No_INBouts = length(IntroNoteResults.NoofINs(find(FirstSyll == 'i')));

IntroNoteResults.MeanIN_NotINBouts = mean(IntroNoteResults.NoofINs(find(FirstSyll ~= 'i')));
IntroNoteResults.STDIN_NotINBouts = std(IntroNoteResults.NoofINs(find(FirstSyll ~= 'i')));
IntroNoteResults.No_NotINBouts = length(IntroNoteResults.NoofINs(find(FirstSyll ~= 'i')));

IntroNoteResults.MeanIN_WithinBouts = mean(IntroNoteResults.WithinBoutNoofINs(find(IntroNoteResults.WithinBoutNoofINs(:,1) > 0),1));
IntroNoteResults.STDIN_WithinBouts = std(IntroNoteResults.WithinBoutNoofINs(find(IntroNoteResults.WithinBoutNoofINs(:,1) > 0),1));
IntroNoteResults.No_WithinBouts = length(IntroNoteResults.WithinBoutNoofINs(find(IntroNoteResults.WithinBoutNoofINs(:,1) > 0),1));
IntroNoteResults.BoutDetails = BoutDetails;


% disp(['Proportion of bouts with an IN before the first motif syllable is ', num2str(length(find(NoofINs > 0))*100/length(NoofINs))]);
% IntroNoteResults.PropINBouts = length(find(NoofINs > 0))*100/length(NoofINs);
% IntroNoteResults.NoofINs = NoofINs;
% 
% 
% % Plots
% 
% % No of intro notes and syllables - distribution
% figure;
% set(gcf, 'Color', 'w');
% Edges = 0:1:max(NoofSylls)+1;
% plot(Edges, histc(NoofINs, Edges)*100/sum(histc(NoofINs, Edges)), 'ks-');
% hold on;
% plot(Edges, histc(NoofSylls, Edges)*100/sum(histc(NoofSylls, Edges)), 'rs-');
% set(gca, 'FontSize', 12);
% xlabel('No of intro notes', 'FontSize', 14);
% ylabel('Percentage of bouts (%)', 'FontSize', 14);
% legend('No of intro notes', 'No of syllables');
% 
% % No of intro notes and time to motif
% figure;
% set(gcf, 'Color', 'w');
% plot(NoofINs, TimeToMotif, 'ks');
% hold on;
% plot(NoofSylls, TimeToMotif, 'rs');
% set(gca, 'FontSize', 12);
% xlabel('No of intro notes', 'FontSize', 14);
% ylabel('Time to Motif (sec)', 'FontSize', 14);
% legend('No of intro notes', 'No of syllables');
% 
% % Corr. between # of sylls / # of INs and bout length
% figure;
% set(gcf, 'Color', 'w');
% subplot(2,2,1);
% plot(NoofINs, [BoutDetails.BoutLength], 'k+');
% ylabel('Bout length (sec)', 'FontSize', 14);
% [r, p] = corrcoef(NoofINs, [BoutDetails.BoutLength]);
% 
% TempAxis = axis;
% text(TempAxis(3) - 1, TempAxis(4) - 1, ['r = ', num2str(r(1,2)), ' and p = ', num2str(p(1,2))]);
% subplot(2,2,2);
% plot(NoofSylls, [BoutDetails.BoutLength], 'k+');
% [r, p] = corrcoef(NoofSylls, [BoutDetails.BoutLength]);
% TempAxis = axis;
% text(TempAxis(3) - 1, TempAxis(4) - 1, ['r = ', num2str(r(1,2)), ' and p = ', num2str(p(1,2))]);
% 
% % Acoustic properties 
% figure;
% set(gcf, 'Color', 'w');
% Feats = [];
% for i = 1:length(BoutDetails),
%     for j = 1:length(INs{i}),
%         if (j == 1)
%             if (j == length(INs{i}))
%                 Feats = [Feats; BoutDetails(i).Feats(INs{i}(j),:) (length(INs{i}) - j + 1)];
%             else
%                 Feats = [Feats; BoutDetails(i).Feats(INs{i}(j),:) 1000];
%             end
%         else
%             Feats = [Feats; BoutDetails(i).Feats(INs{i}(j),:) (length(INs{i}) - j + 1)];
%         end
%     end
% end
% 
% Feats(:,1) = Feats(:,1) * 1000;
% 
% plot3(Feats(find(Feats(:,end) == 1),2), Feats(find(Feats(:,end) == 1),3), Feats(find(Feats(:,end) == 1),1), 'r+', 'MarkerSize', 3);
% hold on;
% 
% MeanLastIN = mean(Feats(find(Feats(:,end) == 1),1:3));
% STDLastIN = std(Feats(find(Feats(:,end) == 1),1:3));
% MeanFirstIN = mean(Feats(find(Feats(:,end) == 1000),1:3));
% STDFirstIN = std(Feats(find(Feats(:,end) == 1000),1:3));
% 
% [x, y, z] = ellipsoid(MeanLastIN(2), MeanLastIN(3), MeanLastIN(1), STDLastIN(2), STDLastIN(3), STDLastIN(1)); 
% surf(x, y, z, 'FaceColor', 'none', 'EdgeColor', 'r');
% 
% plot3(Feats(find(Feats(:,end) == 1000),2), Feats(find(Feats(:,end) == 1000),3), Feats(find(Feats(:,end) == 1000),1), 'b+', 'MarkerSize', 3);
%     
% [x, y, z] = ellipsoid(MeanFirstIN(2), MeanFirstIN(3), MeanFirstIN(1), STDFirstIN(2), STDFirstIN(3), STDFirstIN(1));
% surf(x, y, z, 'FaceColor', 'none', 'EdgeColor', 'b');
% 
% set(gca, 'FontSize', 12);
% xlabel('Log amplitude (dB)', 'FontSize', 14);
% ylabel('Entropy', 'FontSize', 14);
% zlabel('Duration (msec)', 'FontSize', 14);
% legend('Last IN', 'First IN');
% 
% % Acoustic properties 
% figure;
% set(gcf, 'Color', 'w');
% Feats = [];
% for i = 1:length(BoutDetails),
%     for j = 1:length(INs{i}),
%         if (j == 1)
%             if (j == length(INs{i}))
%                 Feats = [Feats; (BoutDetails(i).Feats(INs{i}(j),:)./BoutDetails(i).Feats(INs{i}(end)+1,:)) (length(INs{i}) - j + 1)];
%             else
%                 Feats = [Feats; (BoutDetails(i).Feats(INs{i}(j),:)./BoutDetails(i).Feats(INs{i}(end)+1,:)) 1000];
%             end
%         else
%             Feats = [Feats; (BoutDetails(i).Feats(INs{i}(j),:)./BoutDetails(i).Feats(INs{i}(end)+1,:)) (length(INs{i}) - j + 1)];
%         end
%     end
% end
% 
% Feats(:,1) = Feats(:,1) * 1000;
% plot3(Feats(find(Feats(:,end) == 1),2), Feats(find(Feats(:,end) == 1),3), Feats(find(Feats(:,end) == 1),1), 'r+', 'MarkerSize', 3);
% hold on;
% 
% MeanLastIN = mean(Feats(find(Feats(:,end) == 1),1:3));
% STDLastIN = std(Feats(find(Feats(:,end) == 1),1:3));
% MeanFirstIN = mean(Feats(find(Feats(:,end) == 1000),1:3));
% STDFirstIN = std(Feats(find(Feats(:,end) == 1000),1:3));
% 
% [x, y, z] = ellipsoid(MeanLastIN(2), MeanLastIN(3), MeanLastIN(1), STDLastIN(2), STDLastIN(3), STDLastIN(1)); 
% surf(x, y, z, 'FaceColor', 'none', 'EdgeColor', 'r');
% 
% plot3(Feats(find(Feats(:,end) == 1000),2), Feats(find(Feats(:,end) == 1000),3), Feats(find(Feats(:,end) == 1000),1), 'b+', 'MarkerSize', 3);
%     
% [x, y, z] = ellipsoid(MeanFirstIN(2), MeanFirstIN(3), MeanFirstIN(1), STDFirstIN(2), STDFirstIN(3), STDFirstIN(1));
% surf(x, y, z, 'FaceColor', 'none', 'EdgeColor', 'b');
% 
% set(gca, 'FontSize', 12);
% xlabel('Log amplitude (dB)', 'FontSize', 14);
% ylabel('Entropy', 'FontSize', 14);
% zlabel('Duration (msec)', 'FontSize', 14);
% legend('Last IN', 'First IN');

disp('Finished Analysis');
