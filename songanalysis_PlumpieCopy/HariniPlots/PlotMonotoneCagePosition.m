function [] = PlotMonotoneCagePosition(CSVTextFile, PositionLabels, PositionCoords)

% Position map is as follows:

% Cage is 9" by 9" by 9". The buzzer is placed at the bottom of the cage.
% The bottom is divided into 16 equal parts and the buzzer is placed in
% each of those positions. The mic is roughly in the centre on the top of
% the cage

[BirdParameters, Flag] = ProcessSongData(CSVTextFile, 2000);

% Now for each of the syllables just check that the syllable is complete in
% the file by checking that there is > 50ms before and after the start and
% end of the syllable respectively

Padding = 50;
for i = 1:length(BirdParameters(1).MotifLabels),
    SyllIndices = find(char(BirdParameters(1).SyllableData(:,1)) == BirdParameters(1).MotifLabels(i));
    SyllAmp{i} = [];
    for j = 1:length(SyllIndices),
        if ((BirdParameters(1).SyllableData(SyllIndices(j),4) >= Padding) & (BirdParameters(1).SyllableData(SyllIndices(j),5) <= (BirdParameters(1).FileLen(BirdParameters(1).SyllableData(SyllIndices(j), 2)) - Padding)))
            % SyllAmp{i}(j) = BirdParameters(1).SAPFeatsMatrix(SyllIndices(j),2);
            SyllAmp{i}(j) = BirdParameters(1).AmplitudeKao(SyllIndices(j),1);
        else
            SyllAmp{i}(j) = NaN;
        end
    end
    MeanSyllAmp(i) = nanmean(SyllAmp{i});
    SEMSyllAmp(i) = nanstd(SyllAmp{i})/sqrt(length(find(~isnan(SyllAmp{i}))));
end
                
disp('Finished');