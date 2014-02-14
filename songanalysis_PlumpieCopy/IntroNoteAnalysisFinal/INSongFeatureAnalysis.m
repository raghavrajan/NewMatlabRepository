function [SongFeatureAnalysisResults] = INSongFeatureAnalysis(IntroNoteResults)

MinNumber = 4;

% First calculate feature values for bouts that start with 'i' and have an
% 'i' before the motif

if (size(IntroNoteResults.BoutDetails(1).Feats, 2) == 6)
    FeatureCols = [1 2 3 6];
else
    FeatureCols = [1 2 3];
end
AllINFeats = [];
AllINFeatRatios = [];
AllFirstINFeats = [];% For this one if the number of INs is 1, then the only IN is considered the last and not the first IN
AllLastINFeats = [];
AllFirstINFeats2 = []; % For this one if the number of INs is 1, then the only IN is considered the first

BoutIndices = find((IntroNoteResults.FirstSyll == 'i') & (IntroNoteResults.SyllBeforeMotifs == 'i'));
MaxINs = max(IntroNoteResults.NoofINs(BoutIndices));


for i = 1:MaxINs,
    INBoutFeats{i} = [];
    INBoutIntervals{i} = [];
    for j = 1:length(BoutIndices),
        if (IntroNoteResults.NoofINs(BoutIndices(j)) == i)
            INBoutFeats{i} = [INBoutFeats{i}; [IntroNoteResults.BoutDetails(BoutIndices(j)).Feats([IntroNoteResults.INs{BoutIndices(j)} (IntroNoteResults.INs{BoutIndices(j)}(end) + 1)],FeatureCols) (i:-1:0)']];
            TempFeats = IntroNoteResults.BoutDetails(BoutIndices(j)).Feats(IntroNoteResults.INs{BoutIndices(j)},FeatureCols);
            AllLastINFeats = [AllLastINFeats; [TempFeats(end,:) 1 i]];
            AllFirstINFeats2 = [AllFirstINFeats2; [TempFeats(1,:) 1 i]];
            if (size(TempFeats,1) > 1)
                AllFirstINFeats = [AllFirstINFeats; [TempFeats(1,:) 1 i]];
                for k = 1:(size(TempFeats,1)-1),
                    AllINFeatRatios = [AllINFeatRatios; TempFeats(k+1,:)./TempFeats(k,:)];
                end
            end
        end
    end
end

INBoutFeats_FirstIN = [];
INBoutFeats_SecondLastIN = [];
INBoutFeats_LastIN = [];
INBoutFeats_FirstMotifSyll = [];

for i = 1:MaxINs,
    if (~isempty(INBoutFeats{i}))
        AllINFeats = [AllINFeats; INBoutFeats{i}];
        INBoutFeats_FirstMotifSyll = [INBoutFeats_FirstMotifSyll; INBoutFeats{i}(find(INBoutFeats{i}(:,end) == 0),:)];
        INBoutFeats_LastIN = [INBoutFeats_LastIN; INBoutFeats{i}(find(INBoutFeats{i}(:,end) == 1),:)];
        if (i > 1)
            INBoutFeats_FirstIN = [INBoutFeats_FirstIN; INBoutFeats{i}(find(INBoutFeats{i}(:,end) == i),:)];
        end
    end
end

SongFeatureAnalysisResults.INBoutFeats.FirstMotifSyll = INBoutFeats_FirstMotifSyll;
SongFeatureAnalysisResults.INBoutFeats.LastIN = INBoutFeats_LastIN;
SongFeatureAnalysisResults.INBoutFeats.FirstIN = INBoutFeats_FirstIN;

for i = 1:length(INBoutFeats),
    MeanINBoutFeats{i} = [];
    STDINBoutFeats{i} = [];
    SEMINBoutFeats{i} = [];
    for j = 0:i,
        if (isempty(INBoutFeats{i}))
            continue;
        end
        if (length(find(INBoutFeats{i}(:,end) == j)) > MinNumber)
            TempFeats = INBoutFeats{i}(find(INBoutFeats{i}(:,end) == j),:);
            %TempFeats(:,1:end-1) = TempFeats(:,1:end-1)./INBoutFeats{i}(find(INBoutFeats{i}(:,end) == 0),1:end-1);
            MeanINBoutFeats{i}(j+1,:) = mean(TempFeats);
            STDINBoutFeats{i}(j+1,:) = std(TempFeats);
            SEMINBoutFeats{i}(j+1,:) = std(TempFeats)/sqrt(length(find(INBoutFeats{i}(:,end) == j)));
        end
    end
end

% Now calculate feature values for bouts that start with things other than
% i

BoutIndices = find((IntroNoteResults.FirstSyll ~= 'i') & (IntroNoteResults.SyllBeforeMotifs == 'i'));
MaxINs = max(IntroNoteResults.NoofINs(BoutIndices));

for i = 1:MaxINs,
    Non_INBoutFeats{i} = [];
    for j = 1:length(BoutIndices),
        if (IntroNoteResults.NoofINs(BoutIndices(j)) == i)
            Non_INBoutFeats{i} = [Non_INBoutFeats{i}; [IntroNoteResults.BoutDetails(BoutIndices(j)).Feats([IntroNoteResults.INs{BoutIndices(j)} (IntroNoteResults.INs{BoutIndices(j)}(end) + 1)], FeatureCols) (i:-1:0)']];
            TempFeats = IntroNoteResults.BoutDetails(BoutIndices(j)).Feats(IntroNoteResults.INs{BoutIndices(j)},FeatureCols);
            AllLastINFeats = [AllLastINFeats; [TempFeats(end,:) 2 i]];
            AllFirstINFeats2 = [AllFirstINFeats2; [TempFeats(1,:) 2 i]];
            if (size(TempFeats,1) > 1)
                AllFirstINFeats = [AllFirstINFeats; [TempFeats(1,:) 2 i]];
                for k = 1:(size(TempFeats,1)-1),
                    AllINFeatRatios = [AllINFeatRatios; TempFeats(k+1,:)./TempFeats(k,:)];
                end
            end
        end
    end
end

for i = 1:length(Non_INBoutFeats),
    MeanNon_INBoutFeats{i} = [];
    STDNon_INBoutFeats{i} = [];
    SEMNon_INBoutFeats{i} = [];
    for j = 0:i,
        if (isempty(Non_INBoutFeats{i}))
            continue;
        end
        if (length(find(Non_INBoutFeats{i}(:,end) == j)) > MinNumber)
            TempFeats = Non_INBoutFeats{i}(find(Non_INBoutFeats{i}(:,end) == j),:);
            %TempFeats(:,1:end-1) = TempFeats(:,1:end-1)./Non_INBoutFeats{i}(find(Non_INBoutFeats{i}(:,end) == 0),1:end-1);
            MeanNon_INBoutFeats{i}(j+1,:) = mean(TempFeats);
            STDNon_INBoutFeats{i}(j+1,:) = std(TempFeats);
            SEMNon_INBoutFeats{i}(j+1,:) = std(TempFeats)/sqrt(length(find(Non_INBoutFeats{i}(:,end) == j)));
        end
    end
end

Non_INBoutFeats_FirstIN = [];
Non_INBoutFeats_SecondLastIN = [];
Non_INBoutFeats_LastIN = [];
Non_INBoutFeats_FirstMotifSyll = [];

for i = 1:MaxINs,
    if (~isempty(Non_INBoutFeats{i}))
        AllINFeats = [AllINFeats; Non_INBoutFeats{i}];
        Non_INBoutFeats_FirstMotifSyll = [Non_INBoutFeats_FirstMotifSyll; Non_INBoutFeats{i}(find(Non_INBoutFeats{i}(:,end) == 0),:)];
        Non_INBoutFeats_LastIN = [Non_INBoutFeats_LastIN; Non_INBoutFeats{i}(find(Non_INBoutFeats{i}(:,end) == 1),:)];
        if (i > 1)
            Non_INBoutFeats_FirstIN = [Non_INBoutFeats_FirstIN; Non_INBoutFeats{i}(find(Non_INBoutFeats{i}(:,end) == i),:)];

        end
    end
end

SongFeatureAnalysisResults.Non_INBoutFeats.FirstMotifSyll = Non_INBoutFeats_FirstMotifSyll;
SongFeatureAnalysisResults.Non_INBoutFeats.LastIN = Non_INBoutFeats_LastIN;
SongFeatureAnalysisResults.Non_INBoutFeats.FirstIN = Non_INBoutFeats_FirstIN;

% Now calculate feature values for within bouts

BoutIndices = find(IntroNoteResults.WithinBoutNoofINs(:,1) > 0);
MaxINs = max(IntroNoteResults.WithinBoutNoofINs(BoutIndices,1));

for i = 1:MaxINs,
    Within_INBoutFeats{i} = [];
    for j = 1:length(BoutIndices),
        if (IntroNoteResults.WithinBoutNoofINs(BoutIndices(j),1) == i)
            Within_INBoutFeats{i} = [Within_INBoutFeats{i}; [IntroNoteResults.BoutDetails(IntroNoteResults.WithinBoutINBoutIndices(BoutIndices(j))).Feats([IntroNoteResults.WithinBoutINs{BoutIndices(j)} (IntroNoteResults.WithinBoutINs{BoutIndices(j)}(end) + 1)], FeatureCols) (i:-1:0)']];
            TempFeats = IntroNoteResults.BoutDetails(IntroNoteResults.WithinBoutINBoutIndices(BoutIndices(j))).Feats(IntroNoteResults.WithinBoutINs{BoutIndices(j)}, FeatureCols);
            AllLastINFeats = [AllLastINFeats; [TempFeats(end,:) 3 i]];
            AllFirstINFeats2 = [AllFirstINFeats2; [TempFeats(1,:) 3 i]];
            if (size(TempFeats,1) > 1)
                AllFirstINFeats = [AllFirstINFeats; [TempFeats(1,:) 3 i]];
                for k = 1:(size(TempFeats,1)-1),
                    AllINFeatRatios = [AllINFeatRatios; TempFeats(k+1,:)./TempFeats(k,:)];
                end
            end
        end
    end
end

SongFeatureAnalysisResults.AllINFeatRatios = AllINFeatRatios;
SongFeatureAnalysisResults.AllFirstINFeats = AllFirstINFeats;
SongFeatureAnalysisResults.AllFirstINFeats2 = AllFirstINFeats2;
SongFeatureAnalysisResults.AllLastINFeats = AllLastINFeats;

if (~exist('Within_INBoutFeats', 'var'))
    Within_INBoutFeats = [];
end
    
for i = 1:length(Within_INBoutFeats),
    MeanWithin_INBoutFeats{i} = [];
    STDWithin_INBoutFeats{i} = [];
    SEMWithin_INBoutFeats{i} = [];
    for j = 0:i,
        if (isempty(Within_INBoutFeats{i}))
            continue;
        end
        if (length(find(Within_INBoutFeats{i}(:,end) == j)) > MinNumber)
            TempFeats = Within_INBoutFeats{i}(find(Within_INBoutFeats{i}(:,end) == j),:);
            %TempFeats(:,1:end-1) = TempFeats(:,1:end-1)./Within_INBoutFeats{i}(find(Within_INBoutFeats{i}(:,end) == 0),1:end-1);
            MeanWithin_INBoutFeats{i}(j+1,:) = mean(TempFeats);
            STDWithin_INBoutFeats{i}(j+1,:) = std(TempFeats);
            SEMWithin_INBoutFeats{i}(j+1,:) = std(TempFeats)/sqrt(length(find(Within_INBoutFeats{i}(:,end) == j)));
        end
    end
end

Within_INBoutFeats_FirstIN = [];
Within_INBoutFeats_SecondLastIN = [];
Within_INBoutFeats_LastIN = [];
Within_INBoutFeats_FirstMotifSyll = [];

for i = 1:MaxINs,
    if (~isempty(Within_INBoutFeats{i}))
        AllINFeats = [AllINFeats; Within_INBoutFeats{i}];
        Within_INBoutFeats_FirstMotifSyll = [Within_INBoutFeats_FirstMotifSyll; Within_INBoutFeats{i}(find(Within_INBoutFeats{i}(:,end) == 0),:)];
        Within_INBoutFeats_LastIN = [Within_INBoutFeats_LastIN; Within_INBoutFeats{i}(find(Within_INBoutFeats{i}(:,end) == 1),:)];
        if (i > 1)
            Within_INBoutFeats_FirstIN = [Within_INBoutFeats_FirstIN; Within_INBoutFeats{i}(find(Within_INBoutFeats{i}(:,end) == i),:)];

        end
    end
end

% Now combine all INs across all types - IN bouts, Non IN bouts and within
% bouts based on number of INs in the sequence

MaxINs = max([length(INBoutFeats) length(Non_INBoutFeats) length(Within_INBoutFeats)]);

for i = 1:MaxINs,
    MeanAllINFeats{i} = [];
    STDAllINBoutFeats{i} = [];
    SEMAllINBoutFeats{i} = [];
    CVAllINBoutFeats{i} = [];
    for j = 0:i,
        TempFeats = [];
        
        if (i <= length(INBoutFeats))
            if (~isempty(INBoutFeats{i}))
                TempFeats = [TempFeats; INBoutFeats{i}(find(INBoutFeats{i}(:,end) == j),:)];
            end
        end
        
        if (i <= length(Non_INBoutFeats))
            if (~isempty(Non_INBoutFeats{i}))
                TempFeats = [TempFeats; Non_INBoutFeats{i}(find(Non_INBoutFeats{i}(:,end) == j),:)];
            end
        end
        
        if (i <= length(Within_INBoutFeats))
            if (~isempty(Within_INBoutFeats{i}))
                TempFeats = [TempFeats; Within_INBoutFeats{i}(find(Within_INBoutFeats{i}(:,end) == j),:)];
            end
        end
        
        if (size(TempFeats,1) > MinNumber)
            MeanAllINFeats{i}(j+1,:) = mean(TempFeats);
            STDAllINBoutFeats{i}(j+1,:) = std(TempFeats);
            SEMAllINBoutFeats{i}(j+1,:) = std(TempFeats)/sqrt(size(TempFeats,1));
            CVAllINBoutFeats{i}(j+1,:) = [std(TempFeats(:,1:end-1))./abs(mean(TempFeats(:,1:end-1))) mean(TempFeats(:,end))];
        end
        if (j == 0)
            disp(['Number of INs is ', num2str(i), ' and the number of bouts is ', num2str(size(TempFeats,1))]);
        end
    end
end

SongFeatureAnalysisResults.Within_INBoutFeats.FirstMotifSyll = Within_INBoutFeats_FirstMotifSyll;
SongFeatureAnalysisResults.Within_INBoutFeats.LastIN = Within_INBoutFeats_LastIN;
SongFeatureAnalysisResults.Within_INBoutFeats.FirstIN = Within_INBoutFeats_FirstIN;

SongFeatureAnalysisResults.AllINFeats = AllINFeats;


MatchFlag = 1;
for i = 1:min([length(IntroNoteResults.BoutDetails(1).SongFile) length(IntroNoteResults.BoutDetails(end).SongFile)]);
    MatchFlag = strncmp(IntroNoteResults.BoutDetails(1).SongFile, IntroNoteResults.BoutDetails(end).SongFile, i);
    if (MatchFlag == 0)
        break;
    end
end
TitleString = IntroNoteResults.BoutDetails(1).SongFile(1:i-1);

% Now make the plots
Features = [{'Duration'} {'Amplitude'} {'Entropy'} {'FF'}];
AxisLabels = [{'Full width at half-maximum (sec)'} {'Mean log amplitude (dB)'} {'Wiener Entropy'} {'Fundamental frequency (Hz)'}];
% for Feats = 1:length(Features),
for Feats = 1:length(Features),
    
    figure;
    set(gcf, 'Color', 'w', 'Position', [100 500 350 275]);
    
    %annotation('textbox', [0.05 0.95 0.9 0.025], 'String', [TitleString ': ' Features{Feats}], 'HorizontalAlignment', 'center', 'EdgeColor', 'none')
    
    %subplot(1,2,1);
    hold on;
    for i = 1:length(MeanAllINFeats),
        if (~isempty(MeanAllINFeats{i}))
            if (Feats <= (size(MeanAllINFeats{i}, 2) - 1))
                errorbar(-MeanAllINFeats{i}(2:end,end)+rand(size(MeanAllINFeats{i}(2:end,end)))/5, MeanAllINFeats{i}(2:end,Feats), SEMAllINBoutFeats{i}(2:end,Feats), 'ks-');
            end
        end
    end
%     for i = 1:length(MeanNon_INBoutFeats),
%         if (~isempty(MeanNon_INBoutFeats{i}))
%             errorbar(-MeanNon_INBoutFeats{i}(2:end,end), MeanNon_INBoutFeats{i}(2:end,Feats), SEMNon_INBoutFeats{i}(2:end,Feats), 'rd--');
%         end
%     end
%     for i = 1:length(MeanWithin_INBoutFeats),
%         if (~isempty(MeanWithin_INBoutFeats{i}))
%             errorbar(-MeanWithin_INBoutFeats{i}(2:end,end), MeanWithin_INBoutFeats{i}(2:end,Feats), SEMWithin_INBoutFeats{i}(2:end,Feats), 'bo--');
%         end
%     end
    
%     subplot(1,2,2);
%     hold on;
%     for i = 1:length(MeanINBoutFeats),
%         if (~isempty(MeanINBoutFeats{i}))
%             errorbar(i, MeanINBoutFeats{i}(1,Feats), SEMINBoutFeats{i}(1,Feats), 'ks-');
%         end
%     end
%     for i = 1:length(MeanNon_INBoutFeats),
%         if (~isempty(MeanNon_INBoutFeats{i}))
%             errorbar(i, MeanNon_INBoutFeats{i}(1,Feats), SEMNon_INBoutFeats{i}(1,Feats), 'rd--');
%         end
%     end
%     for i = 1:length(MeanWithin_INBoutFeats),
%         if (~isempty(MeanWithin_INBoutFeats{i}))
%             errorbar(i, MeanWithin_INBoutFeats{i}(1,Feats), SEMWithin_INBoutFeats{i}(1,Feats), 'bo--');
%         end
%     end
    
    set(gca, 'FontSize', 10, 'FontName', 'Arial');
    xlabel('IN position (-1 is closest to motif)', 'FontSize', 12, 'FontName', 'Arial')
    ylabel(AxisLabels{Feats}, 'FontSize', 12, 'FontName', 'Arial');
    axis tight;
    Temp = axis;
    if (Temp(3) > 0)
        NewTemp = [(Temp(1) - 0.2) (Temp(2) + 0.2) (0.99*Temp(3)) (1.01*Temp(4))];
    else
        NewTemp = [(Temp(1) - 0.2) (Temp(2) + 0.2) (1.01*Temp(3)) (0.99*Temp(4))];
    end
    axis(NewTemp);
    
%     saveas(gcf, [TitleString, '.', Features{Feats}, '.fig'], 'fig');
%     saveas(gcf, [TitleString, '.', Features{Feats}, '.png'], 'png');
end

for Feats = 1:length(Features),
    
    figure;
    set(gcf, 'Color', 'w', 'Position', [100 500 350 275]);
    
    %annotation('textbox', [0.05 0.95 0.9 0.025], 'String', [TitleString ': ' Features{Feats}], 'HorizontalAlignment', 'center', 'EdgeColor', 'none')
    
    %subplot(1,2,1);
    hold on;
    for i = 1:length(CVAllINBoutFeats),
        if (~isempty(CVAllINBoutFeats{i}))
            if (Feats <= (size(CVAllINBoutFeats{i}, 2) - 1))
                plot(-CVAllINBoutFeats{i}(2:end,end)+rand(size(MeanAllINFeats{i}(2:end,end)))/5, CVAllINBoutFeats{i}(2:end,Feats), 'ks-');
            end
        end
    end
    
    set(gca, 'FontSize', 10, 'FontName', 'Arial');
    xlabel('IN position (-1 is closest to motif)', 'FontSize', 12, 'FontName', 'Arial')
    ylabel('Co-efficient of variation', 'FontSize', 12, 'FontName', 'Arial');
    axis tight;
    Temp = axis;
    if (Temp(3) > 0)
        NewTemp = [(Temp(1) - 0.2) (Temp(2) + 0.2) (0.99*Temp(3)) (1.01*Temp(4))];
    else
        NewTemp = [(Temp(1) - 0.2) (Temp(2) + 0.2) (1.01*Temp(3)) (0.99*Temp(4))];
    end
    axis(NewTemp);
end

% figure;
% % This figure will have all the trajectory plots for each trial
% subplot(2,2,1); % Amplitude vs. Entropy
% for i = 1:length(INBoutFeats),
%     if (~isempty(INBoutFeats))
%         for j = 1:i+1:size(INBoutFeats{i},1),
%         plot(INBoutFeats{i}(j:(j+i),2), INBoutFeats{i}(j:(j+i),3), 'Color', [0.6 0.6 0.6]);
%         hold on;
%         plot(INBoutFeats{i}(j,2), INBoutFeats{i}(j,3), 'rs');
%         plot(INBoutFeats{i}(j+i,2), INBoutFeats{i}(j+i,3), 'bs');
%         end
%     end
% end
% 
% for i = 1:length(Non_INBoutFeats),
%     if (~isempty(Non_INBoutFeats{i}))
%         for j = 1:i+1:size(Non_INBoutFeats{i},1),
%         plot(Non_INBoutFeats{i}(j:(j+i),2), Non_INBoutFeats{i}(j:(j+i),3), 'c');
%         plot(Non_INBoutFeats{i}(j,2), Non_INBoutFeats{i}(j,3), 'bd', 'MarkerSize', 8);
%         plot(Non_INBoutFeats{i}(j+i,2), Non_INBoutFeats{i}(j+i,3), 'rd', 'MarkerSize', 8);
%         end
%     end
% end
% 
% subplot(2,2,2); % Amplitude vs. Duration
% for i = 1:length(INBoutFeats),
%     if (~isempty(INBoutFeats))
%         for j = 1:i+1:size(INBoutFeats{i},1),
%         plot(INBoutFeats{i}(j:(j+i),2), INBoutFeats{i}(j:(j+i),1), 'Color', [0.6 0.6 0.6]);
%         hold on;
%         plot(INBoutFeats{i}(j,2), INBoutFeats{i}(j,1), 'rs');
%         plot(INBoutFeats{i}(j+i,2), INBoutFeats{i}(j+i,1), 'bs');
%         end
%     end
% end
% 
% for i = 1:length(Non_INBoutFeats),
%     if (~isempty(Non_INBoutFeats{i}))
%         for j = 1:i+1:size(Non_INBoutFeats{i},1),
%         plot(Non_INBoutFeats{i}(j:(j+i),2), Non_INBoutFeats{i}(j:(j+i),1), 'c');
%         plot(Non_INBoutFeats{i}(j,2), Non_INBoutFeats{i}(j,1), 'bd', 'MarkerSize', 8);
%         plot(Non_INBoutFeats{i}(j+i,2), Non_INBoutFeats{i}(j+i,1), 'rd', 'MarkerSize', 8);
%         end
%     end
% end
% 
% subplot(2,2,3); % Entropy vs. Duration
% for i = 1:length(INBoutFeats),
%     if (~isempty(INBoutFeats))
%         for j = 1:i+1:size(INBoutFeats{i},1),
%         plot(INBoutFeats{i}(j:(j+i),3), INBoutFeats{i}(j:(j+i),1), 'Color', [0.6 0.6 0.6]);
%         hold on;
%         plot(INBoutFeats{i}(j,3), INBoutFeats{i}(j,1), 'rs');
%         plot(INBoutFeats{i}(j+i,3), INBoutFeats{i}(j+i,1), 'bs');
%         end
%     end
% end
% 
% for i = 1:length(Non_INBoutFeats),
%     if (~isempty(Non_INBoutFeats{i}))
%         for j = 1:i+1:size(Non_INBoutFeats{i},1),
%         plot(Non_INBoutFeats{i}(j:(j+i),3), Non_INBoutFeats{i}(j:(j+i),1), 'c');
%         plot(Non_INBoutFeats{i}(j,3), Non_INBoutFeats{i}(j,1), 'bd', 'MarkerSize', 8);
%         plot(Non_INBoutFeats{i}(j+i,3), Non_INBoutFeats{i}(j+i,1), 'rd', 'MarkerSize', 8);
%         end
%     end
% end

% figure;
% % all three features
% FirstINBoutFeats = [];
% SecondLastINBoutFeats = [];
% LastINBoutFeats = [];
% MotifFirstSyllBoutFeats = [];
% hold on;
% for i = 1:length(INBoutFeats),
%     if (~isempty(INBoutFeats{i}))
%         if (i > 1)
%             FirstINBoutFeats = [FirstINBoutFeats; INBoutFeats{i}(find(INBoutFeats{i}(:,end) == i),:)];
%         end
%         if (i > 2)
%             SecondLastINBoutFeats = [SecondLastINBoutFeats; INBoutFeats{i}(find(INBoutFeats{i}(:,end) == 2),:)];
%         end
%         LastINBoutFeats = [LastINBoutFeats; INBoutFeats{i}(find(INBoutFeats{i}(:,end) == 1),:)];
%         MotifFirstSyllBoutFeats = [MotifFirstSyllBoutFeats; INBoutFeats{i}(find(INBoutFeats{i}(:,end) == 0),:)];
%         for j = 1:i+1:size(INBoutFeats{i},1),
%             plot3(INBoutFeats{i}(j:(j+i),3), INBoutFeats{i}(j:(j+i),2), INBoutFeats{i}(j:(j+i),1), 'Color', [0.8 0.8 0.8], 'LineWidth', 0.1);
%         end
%     end
% end
% for i = 1:length(INBoutFeats),
%     if (~isempty(INBoutFeats{i}))
%         for j = 1:i+1:size(INBoutFeats{i},1),
%             if (i > 1)
%                 plot3(INBoutFeats{i}(j,3), INBoutFeats{i}(j,2), INBoutFeats{i}(j,1), 'rs', 'MarkerSize', 5);
%             end
%             if (i > 2)
%                 plot3(INBoutFeats{i}(j+i-2,3), INBoutFeats{i}(j+i-2,2), INBoutFeats{i}(j+i-2,1), 'ks', 'MarkerSize', 5);
%             end
%             plot3(INBoutFeats{i}(j+i-1,3), INBoutFeats{i}(j+i-1,2), INBoutFeats{i}(j+i-1,1), 'bs', 'MarkerSize', 5);
%             plot3(INBoutFeats{i}(j+i,3), INBoutFeats{i}(j+i,2), INBoutFeats{i}(j+i,1), 'gs', 'MarkerSize', 5);
%         end
%     end
% end
% MeanFirstINFeats = mean(FirstINBoutFeats);
% STDFirstINFeats = std(FirstINBoutFeats);
% [x, y, z] = ellipsoid(MeanFirstINFeats(3), MeanFirstINFeats(2), MeanFirstINFeats(1), STDFirstINFeats(3), STDFirstINFeats(2), STDFirstINFeats(1));
% FirstINEllipse = surface(x, y, z, 'FaceColor', 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.2);
% 
% MeanSecondLastINFeats = mean(SecondLastINBoutFeats);
% STDSecondLastINFeats = std(SecondLastINBoutFeats);
% [x, y, z] = ellipsoid(MeanSecondLastINFeats(3), MeanSecondLastINFeats(2), MeanSecondLastINFeats(1), STDSecondLastINFeats(3), STDSecondLastINFeats(2), STDSecondLastINFeats(1));
% SecondLastINEllipse = surface(x, y, z, 'FaceColor', 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.2);
% 
% MeanLastINFeats = mean(LastINBoutFeats);
% STDLastINFeats = std(LastINBoutFeats);
% [x, y, z] = ellipsoid(MeanLastINFeats(3), MeanLastINFeats(2), MeanLastINFeats(1), STDLastINFeats(3), STDLastINFeats(2), STDLastINFeats(1));
% LastINEllipse = surface(x, y, z, 'FaceColor', 'b', 'EdgeColor', 'none', 'FaceAlpha', 0.2);
% 
% MeanMotifFirstSyllFeats = mean(MotifFirstSyllBoutFeats);
% STDMotifFirstSyllFeats = std(MotifFirstSyllBoutFeats);
% [x, y, z] = ellipsoid(MeanMotifFirstSyllFeats(3), MeanMotifFirstSyllFeats(2), MeanMotifFirstSyllFeats(1), STDMotifFirstSyllFeats(3), STDMotifFirstSyllFeats(2), STDMotifFirstSyllFeats(1));
% MotifFirstSyllEllipse = surface(x, y, z, 'FaceColor', 'g', 'EdgeColor', 'none', 'FaceAlpha', 0.2);
% 
% set(gca, 'XGrid', 'on')
% set(gca, 'YGrid', 'on')
% set(gca, 'ZGrid', 'on')
% axis tight
% 
% figure;
% hold on;
% for i = 1:length(Non_INBoutFeats),
%     if (~isempty(Non_INBoutFeats{i}))
%         for j = 1:i+1:size(Non_INBoutFeats{i},1),
%             plot3(Non_INBoutFeats{i}(j:(j+i),3), Non_INBoutFeats{i}(j:(j+i),2), Non_INBoutFeats{i}(j:(j+i),1), 'Color', [0.8 0.8 0.8], 'LineWidth', 0.1);
%         end
%     end
% end
% for i = 1:length(Non_INBoutFeats),
%     if (~isempty(Non_INBoutFeats{i}))
%         for j = 1:i+1:size(Non_INBoutFeats{i},1),
%             if (i > 1)
%                 plot3(Non_INBoutFeats{i}(j,3), Non_INBoutFeats{i}(j,2), Non_INBoutFeats{i}(j,1), 'rs', 'MarkerSize', 5);
%             end
%             plot3(Non_INBoutFeats{i}(j+i-1,3), Non_INBoutFeats{i}(j+i-1,2), Non_INBoutFeats{i}(j+i-1,1), 'bs', 'MarkerSize', 5);
%             plot3(Non_INBoutFeats{i}(j+i,3), Non_INBoutFeats{i}(j+i,2), Non_INBoutFeats{i}(j+i,1), 'gs', 'MarkerSize', 5);
%         end
%     end
% end
% [x, y, z] = ellipsoid(MeanFirstINFeats(3), MeanFirstINFeats(2), MeanFirstINFeats(1), STDFirstINFeats(3), STDFirstINFeats(2), STDFirstINFeats(1));
% FirstINEllipse = surface(x, y, z, 'FaceColor', 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.2);
% 
% [x, y, z] = ellipsoid(MeanSecondLastINFeats(3), MeanSecondLastINFeats(2), MeanSecondLastINFeats(1), STDSecondLastINFeats(3), STDSecondLastINFeats(2), STDSecondLastINFeats(1));
% SecondLastINEllipse = surface(x, y, z, 'FaceColor', 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.2);
% 
% [x, y, z] = ellipsoid(MeanLastINFeats(3), MeanLastINFeats(2), MeanLastINFeats(1), STDLastINFeats(3), STDLastINFeats(2), STDLastINFeats(1));
% LastINEllipse = surface(x, y, z, 'FaceColor', 'b', 'EdgeColor', 'none', 'FaceAlpha', 0.2);
% 
% [x, y, z] = ellipsoid(MeanMotifFirstSyllFeats(3), MeanMotifFirstSyllFeats(2), MeanMotifFirstSyllFeats(1), STDMotifFirstSyllFeats(3), STDMotifFirstSyllFeats(2), STDMotifFirstSyllFeats(1));
% MotifFirstSyllEllipse = surface(x, y, z, 'FaceColor', 'g', 'EdgeColor', 'none', 'FaceAlpha', 0.2);
% 
% set(gca, 'XGrid', 'on')
% set(gca, 'YGrid', 'on')
% set(gca, 'ZGrid', 'on')
% axis tight
% 
% 
% figure;
% hold on;
% for i = 1:length(Within_INBoutFeats),
%     if (~isempty(Within_INBoutFeats{i}))
%         for j = 1:i+1:size(Within_INBoutFeats{i},1),
%             plot3(Within_INBoutFeats{i}(j:(j+i),3), Within_INBoutFeats{i}(j:(j+i),2), Within_INBoutFeats{i}(j:(j+i),1), 'Color', [0.8 0.8 0.8], 'LineWidth', 0.1);
%         end
%     end
% end
% for i = 1:length(Within_INBoutFeats),
%     if (~isempty(Within_INBoutFeats{i}))
%         for j = 1:i+1:size(Within_INBoutFeats{i},1),
%             if (i > 1)
%                 plot3(Within_INBoutFeats{i}(j,3), Within_INBoutFeats{i}(j,2), Within_INBoutFeats{i}(j,1), 'rs', 'MarkerSize', 5);
%             end
%             plot3(Within_INBoutFeats{i}(j+i-1,3), Within_INBoutFeats{i}(j+i-1,2), Within_INBoutFeats{i}(j+i-1,1), 'bs', 'MarkerSize', 5);
%             plot3(Within_INBoutFeats{i}(j+i,3), Within_INBoutFeats{i}(j+i,2), Within_INBoutFeats{i}(j+i,1), 'gs', 'MarkerSize', 5);
%         end
%     end
% end
% 
% [x, y, z] = ellipsoid(MeanFirstINFeats(3), MeanFirstINFeats(2), MeanFirstINFeats(1), STDFirstINFeats(3), STDFirstINFeats(2), STDFirstINFeats(1));
% FirstINEllipse = surface(x, y, z, 'FaceColor', 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.2);
% 
% [x, y, z] = ellipsoid(MeanSecondLastINFeats(3), MeanSecondLastINFeats(2), MeanSecondLastINFeats(1), STDSecondLastINFeats(3), STDSecondLastINFeats(2), STDSecondLastINFeats(1));
% SecondLastINEllipse = surface(x, y, z, 'FaceColor', 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.2);
% 
% [x, y, z] = ellipsoid(MeanLastINFeats(3), MeanLastINFeats(2), MeanLastINFeats(1), STDLastINFeats(3), STDLastINFeats(2), STDLastINFeats(1));
% LastINEllipse = surface(x, y, z, 'FaceColor', 'b', 'EdgeColor', 'none', 'FaceAlpha', 0.2);
% 
% [x, y, z] = ellipsoid(MeanMotifFirstSyllFeats(3), MeanMotifFirstSyllFeats(2), MeanMotifFirstSyllFeats(1), STDMotifFirstSyllFeats(3), STDMotifFirstSyllFeats(2), STDMotifFirstSyllFeats(1));
% MotifFirstSyllEllipse = surface(x, y, z, 'FaceColor', 'g', 'EdgeColor', 'none', 'FaceAlpha', 0.2);
% 
% set(gca, 'XGrid', 'on')
% set(gca, 'YGrid', 'on')
% set(gca, 'ZGrid', 'on')
% axis tight


disp('Finished feature analysis');
