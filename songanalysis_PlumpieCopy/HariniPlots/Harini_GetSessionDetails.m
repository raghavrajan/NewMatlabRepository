function [] = Harini_GetSessionDetails(BirdDetailsTextFile)

% First get details from the CSV text file
disp('Getting header data from CSV file ...');
[HeaderLine, BirdDetails] = LSINA_GetDetailsFromCSVFile(BirdDetailsTextFile);

% Now parse all the lines into the appropriate variables based on the
% header line
disp('Getting data from CSV file ...');
[BirdParameters] = LSINA_ParseHeaderBirdData(HeaderLine, BirdDetails);

% Get unique BirdNames to get an idea of how many birds there are
for i = 1:length(BirdParameters),
    BirdNames{i} = BirdParameters(i).BirdName;
    MicrophoneTypes{i} = BirdParameters(i).Microphone;
end

UniqueBirds = unique(BirdNames);

figure;
set(gcf, 'Color', 'w');
set(gcf, 'Position', [318 66 1500 900]);
p = panel();
NumRows = ceil(length(UniqueBirds)/4) + 2;
p.pack(NumRows, 4);
p.fontsize = 10;

Index = 0;
for i = 1:length(UniqueBirds),
    Matches = strmatch(UniqueBirds{i}, BirdNames);
    
    UniqueMicrophones = unique(MicrophoneTypes(Matches));
    
    for MicType = 1:length(UniqueMicrophones),
        Matches = intersect(strmatch(UniqueBirds{i}, BirdNames), strmatch(UniqueMicrophones{MicType}, MicrophoneTypes));
        Index = Index + 1;
        for j = 1:length(Matches),
            Date{Index}{j} = BirdParameters(Matches(j)).DataLabel;
            DateNumber{Index}(j) = datenum(Date{Index}{j}, 'ddmmyy');
            Condition{Index}{j} = BirdParameters(Matches(j)).Condition;
            Microphone{Index}{j} = BirdParameters(Matches(j)).Microphone;
        end
        [SortedVals, SortedIndices] = sort(DateNumber{Index});
        Date{Index} = Date{Index}(SortedIndices);
        DateNumber{Index} = DateNumber{Index}(SortedIndices);
        Condition{Index} = Condition{Index}(SortedIndices);
        Microphone{Index} = Microphone{Index}(SortedIndices);

        p(mod(Index-1,NumRows) + 1, ceil(Index/NumRows)).select();
        for j = 1:length(DateNumber{Index}),
            YVal = [0.8 1.2 1.2 0.8];
            switch Condition{Index}{j}
                case 'L0'
                    YVal = YVal + 0;
                case 'L1'
                    YVal = YVal + 1;
                case 'L2'
                    YVal = YVal + 2;
                case 'L3'
                    YVal = YVal + 3;
                case 'L4'
                    YVal = YVal + 4;
                case 'UN'
                    YVal = YVal + 5;
            end
            switch Microphone{Index}{j}
                case 'MM'
                    PatchBoxColour = 'k';
                case 'HF'
                    PatchBoxColour = 'r';
                case 'BP'
                    PatchBoxColour = 'b';
            end

            patch([(DateNumber{Index}(j)-0.35) (DateNumber{Index}(j)-0.35) (DateNumber{Index}(j)+0.35) (DateNumber{Index}(j)+0.35)] - DateNumber{Index}(1), YVal, 'k', 'EdgeColor', PatchBoxColour, 'FaceColor', PatchBoxColour, 'LineWidth', 2);
        end
        title([UniqueBirds{i} ':' UniqueMicrophones{MicType}]);
        set(gca, 'YTick', 1:1:6, 'YTickLabel', {'L0' 'L1' 'L2' 'L3' 'L4' 'UN'}, 'FontSize', 10);
        
        if (ceil(Index/NumRows) == 1)
            ylabel('Distance');
        end
        if (mod(Index,NumRows) == 0)
            xlabel('Days relative to first day of recording');
        end
        axis tight;
        Temp = axis;
        Temp = [-1 1.05*Temp(2) 0 7];
        axis(Temp);
    end
end
OutputDir = '/home/raghav/StudentRelated/Harini';
set(gcf, 'PaperPositionMode', 'auto');
print(fullfile(OutputDir, 'SessionDatesTimesDetails_SeparateMicrophones.png'), '-dpng');

disp('Finished');