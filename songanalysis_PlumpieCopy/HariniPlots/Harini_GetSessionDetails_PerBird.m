function [] = Harini_GetSessionDetails_PerBird(BirdDetailsTextFile)

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
NumRows = ceil(length(UniqueBirds)/4);
p.pack(NumRows, 4);
p.fontsize = 10;

Index = 1;
for i = 1:length(UniqueBirds),
    Matches = strmatch(UniqueBirds{i}, BirdNames);
    
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

    p(mod(i-1,NumRows) + 1, ceil(i/NumRows)).select();
    for j = 1:length(DateNumber{i}),
        YVal = [0.8 1.2 1.2 0.8];
        switch Condition{i}{j}
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
        switch Microphone{i}{j}
            case 'MM'
                PatchBoxColour = 'k';
            case 'HF'
                PatchBoxColour = 'r';
            case 'BP'
                PatchBoxColour = 'b';
        end

        patch([(DateNumber{i}(j)-0.35) (DateNumber{i}(j)-0.35) (DateNumber{i}(j)+0.35) (DateNumber{i}(j)+0.35)] - DateNumber{i}(1), YVal, 'k', 'EdgeColor', PatchBoxColour, 'FaceColor', PatchBoxColour, 'LineWidth', 2);
    end
    title([UniqueBirds{i}]);
    set(gca, 'YTick', 1:1:6, 'YTickLabel', {'L0' 'L1' 'L2' 'L3' 'L4' 'UN'}, 'FontSize', 10);    
    if (ceil(i/NumRows) == 1)
        ylabel('Distance');
    end
    if (mod(i,NumRows) == 0)
        xlabel('Days relative to first day of recording');
    end
    axis tight;
    Temp = axis;
    Temp = [-1 round(1+Temp(2)) 0 7];
    axis(Temp);
%    set(gca, 'XTick', 0:3:length(Date{i})-1, 'XTickLabel', Date{i}(1:3:end));
%    set(gca, 'XTickLabelRotation', 30);
end
OutputDir = '/home/raghav/StudentRelated/Harini';
set(gcf, 'PaperPositionMode', 'auto');
print(fullfile(OutputDir, 'SessionDatesTimesDetails.png'), '-dpng');
        
disp('Finished');