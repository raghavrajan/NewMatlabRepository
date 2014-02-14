function [PC] = IntroNotePopulationNeuralActivityAnalysisUnWarped(Neural_INR, BinSize, PC)

PreTime = 0.1;
PostTime = 0.2;

Edges = -PreTime:BinSize:PostTime;

IFRFs = 2000;

PreTimeIndex = round(PreTime*IFRFs);
PostTimeIndex = round(PostTime*IFRFs);

Width = 0.005;
GaussianLen = 2;
XGauss = 1:1:(1 + round(2 * GaussianLen * Width * (IFRFs)));
XGauss = XGauss - (length(XGauss) + 1)/2;
GaussWin = (1/((Width * IFRFs) * sqrt(2 * pi))) * exp(-(XGauss.*XGauss)/(2 * (Width * IFRFs) * (Width * IFRFs)));

OverallINIndex = zeros(1, 15);

for NeuronNo = 1:length(Neural_INR),
   
    MaxINs = max([max(Neural_INR{NeuronNo}.NoofINs) max(Neural_INR{NeuronNo}.WithinBoutNoofINs(:,1))]);

    INIndex = zeros(1,MaxINs);

    for i = 1:length(Neural_INR{NeuronNo}.NoofINs),
        if (Neural_INR{NeuronNo}.NoofINs(i) > 0)
            INs = Neural_INR{NeuronNo}.INs{i};
            SpikeTimes = Neural_INR{NeuronNo}.BoutDetails(i).SpikeTimes;
            for j = 1:length(INs),
                INIndex(length(INs) - j + 1) = INIndex(length(INs) - j + 1) + 1;
                Index = INIndex(length(INs) - j + 1);
                
                OverallINIndex(length(INs) - j + 1) = OverallINIndex(length(INs) - j + 1) + 1;
                OverallIndex = OverallINIndex(length(INs) - j + 1);
                
                INOnset = Neural_INR{NeuronNo}.BoutDetails(i).onsets(INs(j));
                INOffset = Neural_INR{NeuronNo}.BoutDetails(i).offsets(INs(j));
                GapOnset = INOffset;
                GapOffset = Neural_INR{NeuronNo}.BoutDetails(i).onsets(INs(j) + 1);

                INDur{NeuronNo}{length(INs) - j + 1}(Index) = INOffset - INOnset;
                GapDur{NeuronNo}{length(INs) - j + 1}(Index) = GapOffset - GapOnset;
    
                AllINDur{length(INs) - j + 1}(OverallIndex) = INOffset - INOnset;
                AllGapDur{length(INs) - j + 1}(OverallIndex) = GapOffset - GapOnset;
                
                TrialSpikeTimes = [];
                TrialSpikeTimes = SpikeTimes(find((SpikeTimes >= (INOnset - PreTime)) & (SpikeTimes < (INOnset + PostTime))));
                PST{NeuronNo}{length(INs) - j + 1}(Index,:) = histc(TrialSpikeTimes, Edges + INOnset)/BinSize;

    %             if (j == length(INs))
    %                 FirstSyllOnset = GapOffset;
    %                 FirstSyllOffset = Neural_INR{NeuronNo}.BoutDetails(i).offsets(INs(j) + 1);
    %                 
    %                 FirstSyllDur(Index) = FirstSyllOffset - FirstSyllOnset;
    %             
    %                 FirstSyllIFROnset = find(IFR(1,:) >= FirstSyllOnset, 1, 'first');
    %                 FirstSyllIFROffset = find(IFR(1,:) >= FirstSyllOffset, 1, 'first');
    %                 
    %                 TempFirstSyllIFRTime = IFR(1,(FirstSyllIFROnset - PreTimeIndex):(FirstSyllIFROffset + PostTimeIndex)) - FirstSyllOnset;
    %                 TempFirstSyllIFR = ST{i}(1, (FirstSyllIFROnset - PreTimeIndex):(FirstSyllIFROffset + PostTimeIndex));
    %                 
    %                 FirstSyllIFR{Index} = [TempFirstSyllIFRTime; TempFirstSyllIFR];
    %             end
            end
        end
    end

    for i = 1:size(Neural_INR{NeuronNo}.WithinBoutNoofINs, 1),
        if (Neural_INR{NeuronNo}.WithinBoutNoofINs(i, 1) > 0)
            INs = Neural_INR{NeuronNo}.WithinBoutINs{i};
            SpikeTimes = Neural_INR{NeuronNo}.BoutDetails(Neural_INR{NeuronNo}.WithinBoutINBoutIndices(i)).SpikeTimes;
            for j = 1:length(INs),
                INIndex(length(INs) - j + 1) = INIndex(length(INs) - j + 1) + 1;
                Index = INIndex(length(INs) - j + 1);

                OverallINIndex(length(INs) - j + 1) = OverallINIndex(length(INs) - j + 1) + 1;
                OverallIndex = OverallINIndex(length(INs) - j + 1);
                
                INOnset = Neural_INR{NeuronNo}.BoutDetails(Neural_INR{NeuronNo}.WithinBoutINBoutIndices(i)).onsets(INs(j));
                INOffset = Neural_INR{NeuronNo}.BoutDetails(Neural_INR{NeuronNo}.WithinBoutINBoutIndices(i)).offsets(INs(j));
                GapOnset = INOffset;
                GapOffset = Neural_INR{NeuronNo}.BoutDetails(Neural_INR{NeuronNo}.WithinBoutINBoutIndices(i)).onsets(INs(j) + 1);
                
                INDur{NeuronNo}{length(INs) - j + 1}(Index) = INOffset - INOnset;
                GapDur{NeuronNo}{length(INs) - j + 1}(Index) = GapOffset - GapOnset;
                
                AllINDur{length(INs) - j + 1}(OverallIndex) = INOffset - INOnset;
                AllGapDur{length(INs) - j + 1}(OverallIndex) = GapOffset - GapOnset;

                TrialSpikeTimes = [];
                TrialSpikeTimes = SpikeTimes(find((SpikeTimes >= (INOnset - PreTime)) & (SpikeTimes < (INOnset + PostTime))));
                PST{NeuronNo}{length(INs) - j + 1}(Index,:) = histc(TrialSpikeTimes, Edges + INOnset)/BinSize;

    %             if (j == length(INs))
    %                 FirstSyllOnset = GapOffset;
    %                 FirstSyllOffset = Neural_INR{NeuronNo}.BoutDetails(Neural_INR{NeuronNo}.WithinBoutINBoutIndices(i)).offsets(INs(j) + 1);
    %                 
    %                 FirstSyllDur(Index) = FirstSyllOffset - FirstSyllOnset;
    %             
    %                 FirstSyllIFROnset = find(IFR(1,:) >= FirstSyllOnset, 1, 'first');
    %                 FirstSyllIFROffset = find(IFR(1,:) >= FirstSyllOffset, 1, 'first');
    %                 
    %                 TempFirstSyllIFRTime = IFR(1,(FirstSyllIFROnset - PreTimeIndex):(FirstSyllIFROffset + PostTimeIndex)) - FirstSyllOnset;
    %                 TempFirstSyllIFR = ST{Neural_INR{NeuronNo}.WithinBoutINBoutIndices(i)}(1, (FirstSyllIFROnset - PreTimeIndex):(FirstSyllIFROffset + PostTimeIndex));
    %                 
    %                 FirstSyllIFR{Index} = [TempFirstSyllIFRTime; TempFirstSyllIFR];
    %             end
            end
        end
    end
end

for i = 1:length(PST),
    for j = 1:length(PST{i}),
        PST{i}{j}(:,end) = [];
    end
end

Edges = Edges(1:end-1) + BinSize/2;

for i = 1:length(PST),
    MeanLastINPST(i,:) = mean(PST{i}{1});
end

if (isempty(PC))
    [PC, Var] = eig(cov(MeanLastINPST'));
    %[PC, Var] = factoran(MeanLastINPST', 3);
end

LastINScore = MeanLastINPST'*PC;

FigNo = get(0, 'CurrentFigure');

if (isempty(FigNo))
    FigNo = 0;
end

figure(FigNo + 1);
hold on;

figure(FigNo + 2);
hold on;

% for i = 1:length(PST),
%     MeanFourthLastINPST(i,:) = mean(PST{i}{4})/std(PST{i}{4}(1:size(PST{i}{4}, 1)*size(PST{i}{4}, 2)));
% end
% 
% FourthLastINScore = MeanFourthLastINPST'*PC;
% 
% plot3(FourthLastINScore(:, end), FourthLastINScore(:, end-1) , FourthLastINScore(:, end-2), 'co-', 'LineWidth', 2);


LastINColor = [1 0.75 0.75];
SecondLastINColor = [0.75 0.75 1];
ThirdLastINColor = [0.75 0.75 0.75];
for j = 1:25,
    clear TempPST TempINDur TempGapDur;
    
    for i = 1:length(PST),
        TrialIndex = randperm(size(PST{i}{1},1));
        TrialIndex = TrialIndex(1);
        TempPST(i,:) = PST{i}{1}(TrialIndex,:);
        TempINDur(i) = INDur{i}{1}(TrialIndex);
        TempGapDur(i) = GapDur{i}{1}(TrialIndex);
    end
    Score = TempPST'*PC;

    INOnsetTime = find(Edges >= 0, 1, 'first');
    GapOnsetTime = find(Edges >= mean(TempINDur), 1, 'first');
    GapOffsetTime = find(Edges >= (mean(TempINDur) + mean(TempGapDur)), 1, 'first');
    if (isempty(GapOffsetTime))
        GapOffsetTime = size(Score,1);
    end
    
    figure(FigNo + 1);
    plot3(Score(1:GapOnsetTime, end), Score(1:GapOnsetTime, end-1) , Score(1:GapOnsetTime, end-2), 'Color', LastINColor);

    plot3(Score(1, end), Score(1, end-1) , Score(1, end-2), 'r+', 'MarkerEdgeColor', LastINColor, 'MarkerFaceColor', LastINColor, 'MarkerSize', 4);
    plot3(Score(INOnsetTime, end), Score(INOnsetTime, end-1) , Score(INOnsetTime, end-2), 'ro', 'MarkerEdgeColor', LastINColor, 'MarkerFaceColor', LastINColor, 'MarkerSize', 4);
    plot3(Score(GapOnsetTime, end), Score(GapOnsetTime, end-1) , Score(GapOnsetTime, end-2), 'r^', 'MarkerEdgeColor', LastINColor, 'MarkerFaceColor', LastINColor, 'MarkerSize', 4);

    figure(FigNo + 2);
    plot3(Score(GapOnsetTime:GapOffsetTime, end), Score(GapOnsetTime:GapOffsetTime, end-1) , Score(GapOnsetTime:GapOffsetTime, end-2), 'Color', LastINColor);

    plot3(Score(1, end), Score(1, end-1) , Score(1, end-2), 'r+', 'MarkerFaceColor', 'r', 'MarkerSize', 4);
    plot3(Score(INOnsetTime, end), Score(INOnsetTime, end-1) , Score(INOnsetTime, end-2), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 4);
    plot3(Score(GapOnsetTime, end), Score(GapOnsetTime, end-1) , Score(GapOnsetTime, end-2), 'r^', 'MarkerFaceColor', 'r', 'MarkerSize', 4);
    
    clear TempPST TempINDur TempGapDur;
    
    for i = 1:length(PST),
        TrialIndex = randperm(size(PST{i}{2},1));
        TrialIndex = TrialIndex(1);
        TempPST(i,:) = PST{i}{2}(TrialIndex,:);
        TempINDur(i) = INDur{i}{2}(TrialIndex);
        TempGapDur(i) = GapDur{i}{2}(TrialIndex);
    end
    Score = TempPST'*PC;

    INOnsetTime = find(Edges >= 0, 1, 'first');
    GapOnsetTime = find(Edges >= mean(TempINDur), 1, 'first');
    GapOffsetTime = find(Edges >= (mean(TempINDur) + mean(TempGapDur)), 1, 'first');
    if (isempty(GapOffsetTime))
        GapOffsetTime = size(Score,1);
    end
    
    figure(FigNo + 1);
    plot3(Score(1:GapOnsetTime, end), Score(1:GapOnsetTime, end-1) , Score(1:GapOnsetTime, end-2), 'Color', SecondLastINColor);

    plot3(Score(1, end), Score(1, end-1) , Score(1, end-2), 'b+', 'MarkerFaceColor', 'b', 'MarkerSize', 4);
    plot3(Score(INOnsetTime, end), Score(INOnsetTime, end-1) , Score(INOnsetTime, end-2), 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 4);
    plot3(Score(GapOnsetTime, end), Score(GapOnsetTime, end-1) , Score(GapOnsetTime, end-2), 'b^', 'MarkerFaceColor', 'b', 'MarkerSize', 4);

    figure(FigNo + 2);
    plot3(Score(GapOnsetTime:GapOffsetTime, end), Score(GapOnsetTime:GapOffsetTime, end-1) , Score(GapOnsetTime:GapOffsetTime, end-2), 'Color', SecondLastINColor);

    plot3(Score(1, end), Score(1, end-1) , Score(1, end-2), 'b+', 'MarkerFaceColor', 'b', 'MarkerSize', 4);
    plot3(Score(INOnsetTime, end), Score(INOnsetTime, end-1) , Score(INOnsetTime, end-2), 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 4);
    plot3(Score(GapOnsetTime, end), Score(GapOnsetTime, end-1) , Score(GapOnsetTime, end-2), 'b^', 'MarkerFaceColor', 'b', 'MarkerSize', 4);
    
    clear TempPST TempINDur TempGapDur;
    
    for i = 1:length(PST),
        TrialIndex = randperm(size(PST{i}{3},1));
        TrialIndex = TrialIndex(1);
        TempPST(i,:) = PST{i}{3}(TrialIndex,:);
        TempINDur(i) = INDur{i}{3}(TrialIndex);
        TempGapDur(i) = GapDur{i}{3}(TrialIndex);
    end
    Score = TempPST'*PC;

    INOnsetTime = find(Edges >= 0, 1, 'first');
    GapOnsetTime = find(Edges >= mean(TempINDur), 1, 'first');
    GapOffsetTime = find(Edges >= (mean(TempINDur) + mean(TempGapDur)), 1, 'first');
    if (isempty(GapOffsetTime))
        GapOffsetTime = size(Score,1);
    end
    
    figure(FigNo + 1);
    plot3(Score(1:GapOnsetTime, end), Score(1:GapOnsetTime, end-1) , Score(1:GapOnsetTime, end-2), 'Color', ThirdLastINColor);

    plot3(Score(1, end), Score(1, end-1) , Score(1, end-2), 'k+', 'MarkerFaceColor', 'k', 'MarkerSize', 4);
    plot3(Score(INOnsetTime, end), Score(INOnsetTime, end-1) , Score(INOnsetTime, end-2), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 4);
    plot3(Score(GapOnsetTime, end), Score(GapOnsetTime, end-1) , Score(GapOnsetTime, end-2), 'k^', 'MarkerFaceColor', 'k', 'MarkerSize', 4);

    figure(FigNo + 2);
    % disp(['Score size:', num2str(size(Score, 1)), ', ', num2str(size(Score, 2)), '; GapOnsetTime:', num2str(GapOnsetTime), '; GapOffsetTime:', num2str(GapOffsetTime)]);
    
    plot3(Score(GapOnsetTime:GapOffsetTime, end), Score(GapOnsetTime:GapOffsetTime, end-1) , Score(GapOnsetTime:GapOffsetTime, end-2), 'Color', ThirdLastINColor);

    plot3(Score(1, end), Score(1, end-1) , Score(1, end-2), 'k+', 'MarkerFaceColor', 'k', 'MarkerSize', 4);
    plot3(Score(INOnsetTime, end), Score(INOnsetTime, end-1) , Score(INOnsetTime, end-2), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 4);
    plot3(Score(GapOnsetTime, end), Score(GapOnsetTime, end-1) , Score(GapOnsetTime, end-2), 'k^', 'MarkerFaceColor', 'k', 'MarkerSize', 4);
    
end

figure(FigNo + 1);
INOnsetTime = find(Edges >= 0, 1, 'first');
GapOnsetTime = find(Edges >= mean(AllINDur{1}), 1, 'first');
GapOffsetTime = find(Edges >= (mean(AllINDur{1}) + mean(AllGapDur{1})), 1, 'first');

plot3(LastINScore(1:GapOnsetTime, end), LastINScore(1:GapOnsetTime, end-1) , LastINScore(1:GapOnsetTime, end-2), 'r-', 'LineWidth', 2);

plot3(LastINScore(1, end), LastINScore(1, end-1) , LastINScore(1, end-2), 'r+', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
plot3(LastINScore(INOnsetTime, end), LastINScore(INOnsetTime, end-1) , LastINScore(INOnsetTime, end-2), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
plot3(LastINScore(GapOnsetTime, end), LastINScore(GapOnsetTime, end-1) , LastINScore(GapOnsetTime, end-2), 'r^', 'MarkerFaceColor', 'r', 'MarkerSize', 8);

figure(FigNo + 2);
plot3(LastINScore(GapOnsetTime:GapOffsetTime, end), LastINScore(GapOnsetTime:GapOffsetTime, end-1) , LastINScore(GapOnsetTime:GapOffsetTime, end-2), 'r-', 'LineWidth', 2);

plot3(LastINScore(1, end), LastINScore(1, end-1) , LastINScore(1, end-2), 'r+', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
plot3(LastINScore(INOnsetTime, end), LastINScore(INOnsetTime, end-1) , LastINScore(INOnsetTime, end-2), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
plot3(LastINScore(GapOnsetTime, end), LastINScore(GapOnsetTime, end-1) , LastINScore(GapOnsetTime, end-2), 'r^', 'MarkerFaceColor', 'r', 'MarkerSize', 8);

for i = 1:length(PST),
    MeanSecondLastINPST(i,:) = mean(PST{i}{2});
end

SecondLastINScore = MeanSecondLastINPST'*PC;

INOnsetTime = find(Edges >= 0, 1, 'first');
GapOnsetTime = find(Edges >= mean(AllINDur{2}), 1, 'first');
GapOffsetTime = find(Edges >= (mean(AllINDur{2}) + mean(AllGapDur{2})), 1, 'first');

figure(FigNo + 1);

plot3(SecondLastINScore(1:GapOnsetTime, end), SecondLastINScore(1:GapOnsetTime, end-1) , SecondLastINScore(1:GapOnsetTime, end-2), 'b-', 'LineWidth', 2);

plot3(SecondLastINScore(1, end), SecondLastINScore(1, end-1) , SecondLastINScore(1, end-2), 'b+', 'MarkerFaceColor', 'b', 'MarkerSize', 8);
plot3(SecondLastINScore(INOnsetTime, end), SecondLastINScore(INOnsetTime, end-1) , SecondLastINScore(INOnsetTime, end-2), 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 8);
plot3(SecondLastINScore(GapOnsetTime, end), SecondLastINScore(GapOnsetTime, end-1) , SecondLastINScore(GapOnsetTime, end-2), 'b^', 'MarkerFaceColor', 'b', 'MarkerSize', 8);

figure(FigNo + 2);
plot3(SecondLastINScore(GapOnsetTime:GapOffsetTime, end), SecondLastINScore(GapOnsetTime:GapOffsetTime, end-1) , SecondLastINScore(GapOnsetTime:GapOffsetTime, end-2), 'b-', 'LineWidth', 2);

plot3(SecondLastINScore(1, end), SecondLastINScore(1, end-1) , SecondLastINScore(1, end-2), 'b+', 'MarkerFaceColor', 'b', 'MarkerSize', 8);
plot3(SecondLastINScore(INOnsetTime, end), SecondLastINScore(INOnsetTime, end-1) , SecondLastINScore(INOnsetTime, end-2), 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 8);
plot3(SecondLastINScore(GapOnsetTime, end), SecondLastINScore(GapOnsetTime, end-1) , SecondLastINScore(GapOnsetTime, end-2), 'b^', 'MarkerFaceColor', 'b', 'MarkerSize', 8);

for i = 1:length(PST),
    MeanThirdLastINPST(i,:) = mean(PST{i}{3});
end

ThirdLastINScore = MeanThirdLastINPST'*PC;

INOnsetTime = find(Edges >= 0, 1, 'first');
GapOnsetTime = find(Edges >= mean(AllINDur{3}), 1, 'first');
GapOffsetTime = find(Edges >= (mean(AllINDur{3}) + mean(AllGapDur{3})), 1, 'first');

figure(FigNo + 1);
plot3(ThirdLastINScore(1:GapOnsetTime, end), ThirdLastINScore(1:GapOnsetTime, end-1) , ThirdLastINScore(1:GapOnsetTime, end-2), 'k', 'LineWidth', 2);

plot3(ThirdLastINScore(1, end), ThirdLastINScore(1, end-1) , ThirdLastINScore(1, end-2), 'k+', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
plot3(ThirdLastINScore(INOnsetTime, end), ThirdLastINScore(INOnsetTime, end-1) , ThirdLastINScore(INOnsetTime, end-2), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
plot3(ThirdLastINScore(GapOnsetTime, end), ThirdLastINScore(GapOnsetTime, end-1) , ThirdLastINScore(GapOnsetTime, end-2), 'k^', 'MarkerFaceColor', 'k', 'MarkerSize', 8);


figure(FigNo + 2);
plot3(ThirdLastINScore(GapOnsetTime:GapOffsetTime, end), ThirdLastINScore(GapOnsetTime:GapOffsetTime, end-1) , ThirdLastINScore(GapOnsetTime:GapOffsetTime, end-2), 'k', 'LineWidth', 2);

plot3(ThirdLastINScore(1, end), ThirdLastINScore(1, end-1) , ThirdLastINScore(1, end-2), 'k+', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
plot3(ThirdLastINScore(INOnsetTime, end), ThirdLastINScore(INOnsetTime, end-1) , ThirdLastINScore(INOnsetTime, end-2), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
plot3(ThirdLastINScore(GapOnsetTime, end), ThirdLastINScore(GapOnsetTime, end-1) , ThirdLastINScore(GapOnsetTime, end-2), 'k^', 'MarkerFaceColor', 'k', 'MarkerSize', 8);

disp('Finished');
