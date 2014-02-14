function [PC] = IntroNotePopulationNeuralActivityAnalysisUnWarped2D(Neural_INR, BinSize, PC)

% Reset random number generator to default state
rng('default');

PreTime = 0.1;
PostTime = 0.2;

Edges = -PreTime:BinSize:PostTime;
BaselineEdges = -0.1:BinSize:-0.05;

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
                %PST{NeuronNo}{length(INs) - j + 1}(Index,:) = histc(TrialSpikeTimes, Edges + INOnset)/BinSize;
                PST{NeuronNo}{length(INs) - j + 1}(Index,:) = histc(TrialSpikeTimes, Edges + INOnset);
                
                if (j == 1)
                    TrialSpikeTimes = SpikeTimes(find((SpikeTimes >= (INOnset + BaselineEdges(1))) & (SpikeTimes < (INOnset + BaselineEdges(end)))));
                    BaselinePST{NeuronNo}(Index,:) = histc(TrialSpikeTimes, BaselineEdges + INOnset)/BinSize;
                end
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
                %PST{NeuronNo}{length(INs) - j + 1}(Index,:) = histc(TrialSpikeTimes, Edges + INOnset)/BinSize;
                PST{NeuronNo}{length(INs) - j + 1}(Index,:) = histc(TrialSpikeTimes, Edges + INOnset);

                if (j == 1)
                    TrialSpikeTimes = SpikeTimes(find((SpikeTimes >= (INOnset + BaselineEdges(1))) & (SpikeTimes < (INOnset + BaselineEdges(end)))));
                    BaselinePST{NeuronNo}(Index,:) = histc(TrialSpikeTimes, BaselineEdges + INOnset)/BinSize;
                end
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

MinNoofINs = min(cellfun(@length, PST));

for i = 1:length(PST),
    MinTrialNos(i) = size(PST{i}{MinNoofINs}, 1);
end

PCAPST = ones(min(MinTrialNos)*MinNoofINs*size(PST{1}{1}, 2), length(PST)) * NaN;

for i = 1:length(PST),
    TempPST = [];
    for j = 1:MinNoofINs,
        for k = 1:min(MinTrialNos),
            TempPST = [TempPST; PST{i}{j}(k,:)'];
        end
    end
    PCAPST(:,i) = TempPST;
end


Edges = Edges(1:end-1) + BinSize/2;

if (isempty(PC))
    [PC, Var] = eig(cov(PCAPST));
    %[PC, Var] = factoran(MeanLastINPST', 3);
end

FigNo = get(0, 'CurrentFigure');

if (isempty(FigNo))
    FigNo = 0;
end

figure(FigNo + 1);
hold on;

figure(FigNo + 2);
hold on;

figure(FigNo + 3);
hold on;

% for i = 1:length(PST),
%     MeanFourthLastINPST(i,:) = mean(PST{i}{4})/std(PST{i}{4}(1:size(PST{i}{4}, 1)*size(PST{i}{4}, 2)));
% end
% 
% FourthLastINScore = MeanFourthLastINPST'*PC;
% 
% plot(FourthLastINScore(:, end), FourthLastINScore(:, end-1) , FourthLastINScore(:, end-2), 'co-', 'LineWidth', 2);


INColors = [1 0.6 0.6; 0.6 0.6 1; 0.6 1 0.6; 0 0.4 0.4; 0.4 0.4 0];
INSymColors = ['rbgcm'];

for INNo = 1:MinNoofINs,
    for j = 1:100,
        clear TempPST TempINDur TempGapDur;
    
        for i = 1:length(PST),
            TrialIndex = randperm(size(PST{i}{INNo},1));
            TrialIndex = TrialIndex(1);
            TempPST(i,:) = PST{i}{INNo}(TrialIndex,:);
            TempINDur(i) = INDur{i}{INNo}(TrialIndex);
            TempGapDur(i) = GapDur{i}{INNo}(TrialIndex);
        end
        Score = TempPST'*PC;

        INScoreTrials{INNo}{j} = Score;

        INOnsetTime = find(Edges >= 0, 1, 'first');
        GapOnsetTime = find(Edges >= mean(TempINDur), 1, 'first');
        GapOffsetTime = find(Edges >= (mean(TempINDur) + mean(TempGapDur)), 1, 'first');
        if (isempty(GapOffsetTime))
            GapOffsetTime = size(Score,1);
        end

%         TrajStartScores{1}(j,:) = Score(1, :);
%         INOnsetTimeScores{1}(j,:) = Score(INOnsetTime, :);
%         GapOnsetTimeScores{1}(j,:) = Score(GapOnsetTime, :);

        figure(FigNo + 1);
        if (j <= 5)
            plot(Score(1:INOnsetTime, end), Score(1:INOnsetTime, end-1), 'Color', INColors(INNo,:));
        end
        plot(Score(1, end), Score(1, end-1) , [INSymColors(INNo), '+'], 'MarkerEdgeColor', INSymColors(INNo), 'MarkerFaceColor', 'w', 'MarkerSize', 2);
        plot(Score(INOnsetTime, end), Score(INOnsetTime, end-1), [INSymColors(INNo), 'o'], 'MarkerEdgeColor', INSymColors(INNo), 'MarkerFaceColor', 'w', 'MarkerSize', 2);
    
        figure(FigNo + 2);
        if (j <= 5)
            plot(Score(INOnsetTime:GapOnsetTime, end), Score(INOnsetTime:GapOnsetTime, end-1), 'Color', INColors(INNo,:));
        end
        plot(Score(INOnsetTime, end), Score(INOnsetTime, end-1), [INSymColors(INNo), 'o'], 'MarkerEdgeColor', INSymColors(INNo), 'MarkerFaceColor', 'w', 'MarkerSize', 2);
        plot(Score(GapOnsetTime, end), Score(GapOnsetTime, end-1), [INSymColors(INNo), '^'], 'MarkerEdgeColor', INSymColors(INNo), 'MarkerFaceColor', 'w', 'MarkerSize', 2);

        figure(FigNo + 3);
        if (j <= 5)
            plot(Score(GapOnsetTime:GapOffsetTime, end), Score(GapOnsetTime:GapOffsetTime, end-1), 'Color', INColors(INNo,:));
        end
        plot(Score(GapOnsetTime, end), Score(GapOnsetTime, end-1), [INSymColors(INNo), '^'], 'MarkerEdgeColor', INSymColors(INNo), 'MarkerFaceColor', 'w', 'MarkerSize', 2);
        plot(Score(GapOffsetTime, end), Score(GapOffsetTime, end-1), [INSymColors(INNo), '+'], 'MarkerEdgeColor', INSymColors(INNo), 'MarkerFaceColor', 'w', 'MarkerSize', 2);    
    end
end

for INNo = 1:MinNoofINs,
    for i = 1:length(PST),
        MeanINPST{INNo}(i,:) = mean(PST{i}{INNo});
    end
    MeanINScore{INNo} = MeanINPST{INNo}'*PC;

    figure(FigNo + 1);

    INOnsetTime = find(Edges >= 0, 1, 'first');
    GapOnsetTime = find(Edges >= mean(AllINDur{INNo}), 1, 'first');
    GapOffsetTime = find(Edges >= (mean(AllINDur{INNo}) + mean(AllGapDur{INNo})), 1, 'first');

    plot(MeanINScore{INNo}(1:INOnsetTime, end), MeanINScore{INNo}(1:INOnsetTime, end-1), INSymColors(INNo), 'LineWidth', 1.5);

    plot(MeanINScore{INNo}(1, end), MeanINScore{INNo}(1, end-1), [INSymColors(INNo), '+'], 'MarkerFaceColor', INSymColors(INNo), 'MarkerSize', 5, 'LineWidth', 2);
    plot(MeanINScore{INNo}(INOnsetTime, end), MeanINScore{INNo}(INOnsetTime, end-1), [INSymColors(INNo), 'o'], 'MarkerFaceColor', INSymColors(INNo), 'MarkerSize', 5);

    figure(FigNo + 2);
    plot(MeanINScore{INNo}(INOnsetTime:GapOnsetTime, end), MeanINScore{INNo}(INOnsetTime:GapOnsetTime, end-1), INSymColors(INNo), 'LineWidth', 1.5);

    plot(MeanINScore{INNo}(INOnsetTime, end), MeanINScore{INNo}(INOnsetTime, end-1), [INSymColors(INNo), 'o'], 'MarkerFaceColor', INSymColors(INNo), 'MarkerSize', 5);
    plot(MeanINScore{INNo}(GapOnsetTime, end), MeanINScore{INNo}(GapOnsetTime, end-1), [INSymColors(INNo), '^'], 'MarkerFaceColor', INSymColors(INNo), 'MarkerSize', 5);
    
    figure(FigNo + 3);
    plot(MeanINScore{INNo}(GapOnsetTime:GapOffsetTime, end), MeanINScore{INNo}(GapOnsetTime:GapOffsetTime, end-1), INSymColors(INNo), 'LineWidth', 1.5);

    plot(MeanINScore{INNo}(GapOnsetTime, end), MeanINScore{INNo}(GapOnsetTime, end-1), [INSymColors(INNo), '^'], 'MarkerFaceColor', INSymColors(INNo), 'MarkerSize', 5);
    plot(MeanINScore{INNo}(GapOffsetTime, end), MeanINScore{INNo}(GapOffsetTime, end-1), [INSymColors(INNo), '+'], 'MarkerFaceColor', INSymColors(INNo), 'MarkerSize', 5, 'LineWidth', 2);
end

diag(Var)/sum(diag(Var)) * 100

for i = 1:MinNoofINs,
    INDist{i} = [];
    for j = 1:100,
        if (i == 1)
            for k = j+1:100,
               INDist{i} = [INDist{i} diag(pdist2(INScoreTrials{1}{j}, INScoreTrials{i}{k}))];
            end
        else
            for k = 1:100,
                INDist{i} = [INDist{i} diag(pdist2(INScoreTrials{1}{j}, INScoreTrials{i}{k}))];
            end
        end
    end
end

figure(FigNo + 4);
hold on;
for i = 1:MinNoofINs,
    fill([Edges fliplr(Edges)], [(mean(INDist{i}') + std(INDist{i}')) fliplr((mean(INDist{i}') - std(INDist{i}')))], 'k', 'EdgeColor', 'none', 'FaceColor', INColors(i,:), 'FaceAlpha', 0.5);
    plot(Edges, mean(INDist{i}'), [INSymColors(i), 'o-'], 'MarkerSize', 3);
end

axis tight;
tempaxis = axis;
tempaxis = [-0.1 0.2 0 tempaxis(4)*1.05];
axis(tempaxis);

for i = 1:MinNoofINs,
    plot([mean(AllINDur{i}) mean(AllINDur{i})], [0 tempaxis(4)], [INSymColors(i), '--']);
    plot([(mean(AllINDur{i}) + mean(AllGapDur{i})) (mean(AllINDur{i}) + mean(AllGapDur{i}))], [0 tempaxis(4)], [INSymColors(i), ':']);
end

for i = 1:MinNoofINs,
    INPosDist{i} = [];
    for j = 1:100,
        for k = j+1:100,
            INPosDist{i} = [INPosDist{i} diag(pdist2(INScoreTrials{i}{j}, INScoreTrials{i}{k}))];
        end
    end
end

figure(FigNo + 5);
hold on;
for i = 1:MinNoofINs,
    fill([Edges fliplr(Edges)], [(mean(INPosDist{i}') + std(INPosDist{i}')) fliplr((mean(INPosDist{i}') - std(INPosDist{i}')))], 'k', 'EdgeColor', 'none', 'FaceColor', INColors(i,:), 'FaceAlpha', 0.5);
    plot(Edges, mean(INPosDist{i}'), [INSymColors(i), 'o-'], 'MarkerSize', 3);
end

axis tight;
tempaxis = axis;
tempaxis = [-0.1 0.2 0 tempaxis(4)*1.05];
axis(tempaxis);

for i = 1:MinNoofINs,
    plot([mean(AllINDur{i}) mean(AllINDur{i})], [0 tempaxis(4)], [INSymColors(i), '--']);
    plot([(mean(AllINDur{i}) + mean(AllGapDur{i})) (mean(AllINDur{i}) + mean(AllGapDur{i}))], [0 tempaxis(4)], [INSymColors(i), ':']);
end
disp('Finished');
