function [] = Harini_SessionWiseAmplitueAnalysis(IndividualBirds, BirdNames, BirdOption, FeatureToPlot_ColNo, PlotType, DifferentMicrophones)

Colours = 'rgbcmk';
Symbols = '+o<sd';


FinalFigureDir = '/home/raghav/StudentRelated/Harini/Manuscript/results';

Colours = distinguishable_colors(6); % for the different contexts
Symbols = 'so^'; % for D, DUN and UN songs
for Bird = 1:length(IndividualBirds),
    Sylls = IndividualBirds(Bird).SortedBirdParameters(1).MotifLabels;
    figure;
    p = panel();
    p.pack(2, ceil(length(Sylls)/2));
        
    Conditions = unique(IndividualBirds(Bird).AllConditionIndices);
    RecordingDays = unique(IndividualBirds(Bird).AllRecordingDayIndices);
    for i = RecordingDays(:)',
        for j = Conditions(:)',
            for Syll = 1:length(Sylls),
                p(mod(Syll-1, 2)+1, ceil(Syll/2)).select();
                hold on;
                Matches = find((char(IndividualBirds(Bird).AllSyllableData(:,1)) == Sylls(Syll)) & (IndividualBirds(Bird).AllRecordingDayIndices(:) == i) & (IndividualBirds(Bird).AllConditionIndices(:) == j));
                if (~isempty(Matches))
                    Amplitudes = IndividualBirds(Bird).AdjustedSyllLogAmplitudeMeanValue(Matches);
                    SongType = IndividualBirds(Bird).AllSyllableCategorisation(Matches);
                    SyllTime = IndividualBirds(Bird).AllSyllableTime(Matches);
                    
                    if (j == 6)
                        SongIndices = 1:1:length(Matches);
                        if (~isempty(SongIndices))
                            errorbar((i-1)*24 + mean(SyllTime) + ((2)*10)/60, mean(Amplitudes(SongIndices)), std(Amplitudes(SongIndices))/sqrt(length(SongIndices)), 'ks-', 'Color', Colours(j,:), 'Marker', Symbols(3), 'MarkerSize', 6, 'MarkerFaceColor', Colours(j,:));
                        end
                    else
                        % Now to divide it into the different song types
                        SongTypes = {'D' 'DUN' 'UN'};
                        for k = 1:length(SongTypes),
                            SongIndices = strmatch(SongTypes{k}, SongType, 'exact');
                            if (~isempty(SongIndices))
                                errorbar((i-1)*24 + mean(SyllTime) + ((k-1)*10)/60, mean(Amplitudes(SongIndices)), std(Amplitudes(SongIndices))/sqrt(length(SongIndices)), 'ks-', 'Color', Colours(j,:), 'Marker', Symbols(k), 'MarkerSize', 6, 'MarkerFaceColor', Colours(j,:));
                            end
                        end        
                    end
                    MeanSyllTime{i}(j,Syll) = mean(SyllTime);
                else
                    MeanSyllTime{i}(j,Syll) = NaN;
                end
            end
        end
    end
end
disp('Finished');