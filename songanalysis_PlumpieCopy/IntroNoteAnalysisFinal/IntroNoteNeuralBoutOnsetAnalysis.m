function [] = IntroNoteNeuralBoutOnsetAnalysis(Neural_INR, Colour, FillColour, Colour2, FillColour2)

BinSize = 0.05;
PlotEdges = -1.5:BinSize:0.25;
TransPoint = find(PlotEdges <= -1, 1, 'last');

for i = 1:length(Neural_INR),
	for j = 1:length(Neural_INR{i}.BoutDetails),
        Edges = PlotEdges + Neural_INR{i}.BoutDetails(j).onsets(1);
        PST{i}(j,:) = histc(Neural_INR{i}.BoutDetails(j).SpikeTimes, Edges);
    end
    MeanPST(i,:) = mean(PST{i})/BinSize;
    MeanNormPST(i,:) = mean(PST{i})/BinSize - mean(mean(PST{i}(:,1:TransPoint))/BinSize);
end
MeanPST(:,end) = [];
MeanNormPST(:,end) = [];

Index = 0;
for i = 1:length(Neural_INR),
    if ((length(Neural_INR{i}.AllBoutDetails) - length(Neural_INR{i}.BoutDetails)) > 2)
        Index = Index + 1;
        for j = 1:length(Neural_INR{i}.BoutDetails),
            Edges = PlotEdges + Neural_INR{i}.BoutDetails(j).onsets(1);
            MatchedSongPST{Index}(j,:) = histc(Neural_INR{i}.BoutDetails(j).SpikeTimes, Edges);
        end
        MeanMatchedSongPST(Index,:) = mean(MatchedSongPST{Index})/BinSize;
        MeanNormMatchedSongPST(Index,:) = mean(MatchedSongPST{Index})/BinSize - mean(mean(MatchedSongPST{Index}(:,1:TransPoint))/BinSize);
        
        Indices = find([Neural_INR{i}.AllBoutDetails.Motif] == 0);
        RowNo = 0;
        for j = Indices,
            RowNo = RowNo + 1;
            Edges = PlotEdges + Neural_INR{i}.AllBoutDetails(j).onsets(1);
            MatchedNonSongPST{Index}(RowNo,:) = histc(Neural_INR{i}.AllBoutDetails(j).SpikeTimes, Edges);
        end
        MeanMatchedNonSongPST(Index,:) = mean(MatchedNonSongPST{Index})/BinSize;
        MeanNormMatchedNonSongPST(Index,:) = mean(MatchedNonSongPST{Index})/BinSize - mean(mean(MatchedNonSongPST{Index}(:,1:TransPoint))/BinSize);
    end
end
MeanNormMatchedSongPST(:,end) = [];
MeanNormMatchedNonSongPST(:,end) = [];

PlotEdges(end) = [];

figure;
hold on;
fill([PlotEdges fliplr(PlotEdges)], [(mean(MeanNormPST) - std(MeanNormPST)/sqrt(size(MeanNormPST,1))) fliplr((mean(MeanNormPST) + std(MeanNormPST)/sqrt(size(MeanNormPST,1))))], FillColour, 'EdgeColor', FillColour)
plot(PlotEdges, mean(MeanNormPST), Colour);

plot([PlotEdges(1) PlotEdges(end)], [(mean(mean(MeanNormPST(:,1:TransPoint))) + 2*std(mean(MeanNormPST(:,1:TransPoint)))) (mean(mean(MeanNormPST(:,1:TransPoint))) + 2*std(mean(MeanNormPST(:,1:TransPoint))))], 'k--', 'LineWidth', 1);
plot([PlotEdges(1) PlotEdges(end)], [(mean(mean(MeanNormPST(:,1:TransPoint))) - 2*std(mean(MeanNormPST(:,1:TransPoint)))) (mean(mean(MeanNormPST(:,1:TransPoint))) - 2*std(mean(MeanNormPST(:,1:TransPoint))))], 'k--', 'LineWidth', 1);

fill([0 0.1 0.1 0], [-15 -15 250 250], [0.85 0.85 0.85], 'FaceAlpha', 0.5, 'FaceColor', [0.85 0.85 0.85], 'EdgeColor', 'none');

figure;
hold on;
fill([PlotEdges fliplr(PlotEdges)], [(mean(MeanNormMatchedSongPST) - std(MeanNormMatchedSongPST)/sqrt(size(MeanNormMatchedSongPST,1))) fliplr((mean(MeanNormMatchedSongPST) + std(MeanNormMatchedSongPST)/sqrt(size(MeanNormMatchedSongPST,1))))], FillColour, 'EdgeColor', FillColour)
plot(PlotEdges, mean(MeanNormMatchedSongPST), Colour);

fill([PlotEdges fliplr(PlotEdges)], [(mean(MeanNormMatchedNonSongPST) - std(MeanNormMatchedNonSongPST)/sqrt(size(MeanNormMatchedNonSongPST,1))) fliplr((mean(MeanNormMatchedNonSongPST) + std(MeanNormMatchedNonSongPST)/sqrt(size(MeanNormMatchedNonSongPST,1))))], FillColour2, 'EdgeColor', FillColour2)
plot(PlotEdges, mean(MeanNormMatchedNonSongPST), Colour2);

plot([PlotEdges(1) PlotEdges(end)], [(mean(mean(MeanNormMatchedSongPST(:,1:TransPoint))) + 2*std(mean(MeanNormMatchedSongPST(:,1:TransPoint)))) (mean(mean(MeanNormPST(:,1:TransPoint))) + 2*std(mean(MeanNormPST(:,1:TransPoint))))], 'k--', 'LineWidth', 2);
plot([PlotEdges(1) PlotEdges(end)], [(mean(mean(MeanNormMatchedSongPST(:,1:TransPoint))) - 2*std(mean(MeanNormMatchedSongPST(:,1:TransPoint)))) (mean(mean(MeanNormPST(:,1:TransPoint))) - 2*std(mean(MeanNormPST(:,1:TransPoint))))], 'k--', 'LineWidth', 2);

disp('Finished Analysis');