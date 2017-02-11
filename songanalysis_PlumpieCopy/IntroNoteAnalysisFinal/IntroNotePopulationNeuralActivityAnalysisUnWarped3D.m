function [PC] = IntroNotePopulationNeuralActivityAnalysisUnWarped3D(Neural_INR, BinSize, PC)

% Reset random number generator to default state
rng('default');

PreTime = 0.3;
PostTime = 0.3;

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
    
    for i = 1:MaxINs,
        Raster{NeuronNo}{i} = [];
    end
    
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
                Raster{NeuronNo}{length(INs) - j + 1} = [Raster{NeuronNo}{length(INs) - j + 1}; [(TrialSpikeTimes(:) - INOnset) ones(size(TrialSpikeTimes(:)))*Index]];
                
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
                Raster{NeuronNo}{length(INs) - j + 1} = [Raster{NeuronNo}{length(INs) - j + 1}; [(TrialSpikeTimes(:) - INOnset) ones(size(TrialSpikeTimes(:)))*Index]];
                
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
    TempPST = [];
    for j = 1:length(PST{i}),
        PST{i}{j}(:,end) = [];
        TempPST = [TempPST; PST{i}{j}(:)];
    end
    for j = 1:length(PST{i}),
        PST{i}{j} = (PST{i}{j} - mean(TempPST))/std(TempPST);
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


INColors = [1 0.75 0.75; 0.75 0.75 1; 0.75 1 0.75; 0 0.4 0.4; 0.4 0.4 0];
INSymColors = ['rbgcm'];

NoofReps = 200;

NoofTrajectories = 10;
NoofPoints = 50;

for INNo = 1:MinNoofINs,
    for j = 1:NoofReps,
        clear TempPST TempINDur TempGapDur;
    
        for i = 1:length(PST),
            TrialIndex = randperm(size(PST{i}{INNo},1));
            TrialIndex = TrialIndex(1);
            TempPST(i,:) = PST{i}{INNo}(TrialIndex,:);
            TempINDur(i) = INDur{i}{INNo}(TrialIndex);
            TempGapDur(i) = GapDur{i}{INNo}(TrialIndex);
        end
        Score = TempPST'*PC;

        INTrials{INNo}{j} = TempPST';
        INScoreTrials{INNo}{j} = Score;
        INDurTrials{INNo}{j} = TempINDur;
        INGapDurTrials{INNo}{j} = TempGapDur;
    end
end

for INNo = 1:1,
    for j = 1:NoofTrajectories,
        Score = INTrials{INNo}{j}*PC;
        
        INOnsetTime = find(Edges >= 0, 1, 'first');
        GapOnsetTime = find(Edges >= mean(INDurTrials{INNo}{j}), 1, 'first');
        GapOffsetTime = find(Edges >= (mean(INDurTrials{INNo}{j}) + mean(INGapDurTrials{INNo}{j})), 1, 'first');
        if (isempty(GapOffsetTime))
            GapOffsetTime = size(Score,1);
        end

        figure(FigNo + 1);
        plot3(Score(1:INOnsetTime, end), Score(1:INOnsetTime, end-1), Score(1:INOnsetTime, end-2), 'Color', INColors(INNo,:));
 
        figure(FigNo + 2);
        plot3(Score(INOnsetTime:GapOnsetTime, end), Score(INOnsetTime:GapOnsetTime, end-1), Score(INOnsetTime:GapOnsetTime, end-2), 'Color', INColors(INNo,:));
        
        figure(FigNo + 3);
        plot3(Score(GapOnsetTime:GapOffsetTime, end), Score(GapOnsetTime:GapOffsetTime, end-1), Score(GapOnsetTime:GapOffsetTime, end-2), 'Color', INColors(INNo,:));
    end
end

for INNo = 1:1,
    for j = 1:NoofPoints,
        Score = INTrials{INNo}{j}*PC;
        
        INOnsetTime = find(Edges >= 0, 1, 'first');
        GapOnsetTime = find(Edges >= mean(INDurTrials{INNo}{j}), 1, 'first');
        GapOffsetTime = find(Edges >= (mean(INDurTrials{INNo}{j}) + mean(INGapDurTrials{INNo}{j})), 1, 'first');
        if (isempty(GapOffsetTime))
            GapOffsetTime = size(Score,1);
        end
        
        figure(FigNo + 1);
        plot3(Score(1, end), Score(1, end-1), Score(1, end-2), [INSymColors(INNo), '+'], 'MarkerEdgeColor', INSymColors(INNo), 'MarkerFaceColor', 'w', 'MarkerSize', 2);
        plot3(Score(INOnsetTime, end), Score(INOnsetTime, end-1), Score(INOnsetTime, end-2), [INSymColors(INNo), 'o'], 'MarkerEdgeColor', INSymColors(INNo), 'MarkerFaceColor', 'w', 'MarkerSize', 2);

        figure(FigNo + 2);
        plot3(Score(INOnsetTime, end), Score(INOnsetTime, end-1), Score(INOnsetTime, end-2), [INSymColors(INNo), 'o'], 'MarkerEdgeColor', INSymColors(INNo), 'MarkerFaceColor', 'w', 'MarkerSize', 2);
        plot3(Score(GapOnsetTime, end), Score(GapOnsetTime, end-1), Score(GapOnsetTime, end-2), [INSymColors(INNo), '^'], 'MarkerEdgeColor', INSymColors(INNo), 'MarkerFaceColor', 'w', 'MarkerSize', 2);

        figure(FigNo + 3);
        plot3(Score(GapOnsetTime, end), Score(GapOnsetTime, end-1), Score(GapOnsetTime, end-2), [INSymColors(INNo), '^'], 'MarkerEdgeColor', INSymColors(INNo), 'MarkerFaceColor', 'w', 'MarkerSize', 2);
        plot3(Score(GapOffsetTime, end), Score(GapOffsetTime, end-1), Score(GapOffsetTime, end-2), [INSymColors(INNo), '+'], 'MarkerEdgeColor', INSymColors(INNo), 'MarkerFaceColor', 'w', 'MarkerSize', 2);    
    end
end
    
for INNo = 1:1,
    for i = 1:length(PST),
        MeanINPST{INNo}(i,:) = mean(PST{i}{INNo});
    end
    MeanINScore{INNo} = MeanINPST{INNo}'*PC;

    figure(FigNo + 1);

    INOnsetTime = find(Edges >= 0, 1, 'first');
    GapOnsetTime = find(Edges >= mean(AllINDur{INNo}), 1, 'first');
    GapOffsetTime = find(Edges >= (mean(AllINDur{INNo}) + mean(AllGapDur{INNo})), 1, 'first');

    plot3(MeanINScore{INNo}(1:INOnsetTime, end), MeanINScore{INNo}(1:INOnsetTime, end-1), MeanINScore{INNo}(1:INOnsetTime, end-2), INSymColors(INNo), 'LineWidth', 2);

    plot3(MeanINScore{INNo}(1, end), MeanINScore{INNo}(1, end-1), MeanINScore{INNo}(1, end-2), [INSymColors(INNo), '+'], 'MarkerFaceColor', INSymColors(INNo), 'MarkerSize', 5, 'LineWidth', 3);
    plot3(MeanINScore{INNo}(INOnsetTime, end), MeanINScore{INNo}(INOnsetTime, end-1), MeanINScore{INNo}(INOnsetTime, end-2), [INSymColors(INNo), 'o'], 'MarkerFaceColor', INSymColors(INNo), 'MarkerSize', 5);

    figure(FigNo + 2);
    plot3(MeanINScore{INNo}(INOnsetTime:GapOnsetTime, end), MeanINScore{INNo}(INOnsetTime:GapOnsetTime, end-1), MeanINScore{INNo}(INOnsetTime:GapOnsetTime, end-2), INSymColors(INNo), 'LineWidth', 2);

    plot3(MeanINScore{INNo}(INOnsetTime, end), MeanINScore{INNo}(INOnsetTime, end-1), MeanINScore{INNo}(INOnsetTime, end-2), [INSymColors(INNo), 'o'], 'MarkerFaceColor', INSymColors(INNo), 'MarkerSize', 5);
    plot3(MeanINScore{INNo}(GapOnsetTime, end), MeanINScore{INNo}(GapOnsetTime, end-1), MeanINScore{INNo}(GapOnsetTime, end-2), [INSymColors(INNo), '^'], 'MarkerFaceColor', INSymColors(INNo), 'MarkerSize', 5);
    
    figure(FigNo + 3);
    plot3(MeanINScore{INNo}(GapOnsetTime:GapOffsetTime, end), MeanINScore{INNo}(GapOnsetTime:GapOffsetTime, end-1), MeanINScore{INNo}(GapOnsetTime:GapOffsetTime, end-2), INSymColors(INNo), 'LineWidth', 2);

    plot3(MeanINScore{INNo}(GapOnsetTime, end), MeanINScore{INNo}(GapOnsetTime, end-1), MeanINScore{INNo}(GapOnsetTime, end-2), [INSymColors(INNo), '^'], 'MarkerFaceColor', INSymColors(INNo), 'MarkerSize', 5);
    plot3(MeanINScore{INNo}(GapOffsetTime, end), MeanINScore{INNo}(GapOffsetTime, end-1), MeanINScore{INNo}(GapOffsetTime, end-2), [INSymColors(INNo), '+'], 'MarkerFaceColor', INSymColors(INNo), 'MarkerSize', 5, 'LineWidth', 3);
end

diag(Var)/sum(diag(Var)) * 100

PCComps = [round(size(PC,2)/2):size(PC,2)];

for i = 1:MinNoofINs,
    INDist{i} = [];
    for j = 1:NoofReps,
        if (i == 1)
            for k = j+1:NoofReps,
               INDist{i} = [INDist{i} diag(pdist2(INTrials{1}{j}, INTrials{i}{k}))];
            end
        else
            for k = 1:NoofReps,
                INDist{i} = [INDist{i} diag(pdist2(INTrials{1}{j}, INTrials{i}{k}))];
            end
        end
    end
end

figure(FigNo + 4);
hold on;
for i = 1:MinNoofINs,
    fill([Edges fliplr(Edges)], [(mean(INDist{i}') + 1*std(INDist{i}')) fliplr((mean(INDist{i}') - 1*std(INDist{i}')))], 'k', 'EdgeColor', 'none', 'FaceColor', INColors(i,:), 'FaceAlpha', 0.5);
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
    for j = 1:NoofReps,
        for k = j+1:NoofReps,
            INPosDist{i} = [INPosDist{i} diag(pdist2(INTrials{i}{j}, INTrials{i}{k}))];
        end
    end
end

figure(FigNo + 5);
hold on;
for i = 1:MinNoofINs,
    fill([Edges fliplr(Edges)], [(mean(INPosDist{i}') + std(INPosDist{i}')) fliplr((mean(INPosDist{i}') - std(INPosDist{i}')))], 'k', 'EdgeColor', 'none', 'FaceColor', INColors(i,:), 'FaceAlpha', 0.5);
    plot(Edges, mean(INPosDist{i}'), [INSymColors(i), 'o-'], 'MarkerSize', 3);
end

axis(tempaxis);

for i = 1:MinNoofINs,
    plot([mean(AllINDur{i}) mean(AllINDur{i})], [0 tempaxis(4)], [INSymColors(i), '--']);
    plot([(mean(AllINDur{i}) + mean(AllGapDur{i})) (mean(AllINDur{i}) + mean(AllGapDur{i}))], [0 tempaxis(4)], [INSymColors(i), ':']);
end

for i = 1:MinNoofINs,
    INSecondLastDist{i} = [];
    for j = 1:NoofReps,
        if (i == 2)
            for k = j+1:NoofReps,
               INSecondLastDist{i} = [INSecondLastDist{i} diag(pdist2(INTrials{2}{j}, INTrials{i}{k}))];
            end
        else
            for k = 1:NoofReps,
                INSecondLastDist{i} = [INSecondLastDist{i} diag(pdist2(INTrials{2}{j}, INTrials{i}{k}))];
            end
        end
    end
end

figure(FigNo + 6);
hold on;
for i = 1:MinNoofINs,
    fill([Edges fliplr(Edges)], [(mean(INSecondLastDist{i}') + 2*std(INSecondLastDist{i}')) fliplr((mean(INSecondLastDist{i}') - 2*std(INSecondLastDist{i}')))], 'k', 'EdgeColor', 'none', 'FaceColor', INColors(i,:), 'FaceAlpha', 0.5);
    plot(Edges, mean(INSecondLastDist{i}'), [INSymColors(i), 'o-'], 'MarkerSize', 3);
end

axis(tempaxis);

for i = 1:MinNoofINs,
    plot([mean(AllINDur{i}) mean(AllINDur{i})], [0 tempaxis(4)], [INSymColors(i), '--']);
    plot([(mean(AllINDur{i}) + mean(AllGapDur{i})) (mean(AllINDur{i}) + mean(AllGapDur{i}))], [0 tempaxis(4)], [INSymColors(i), ':']);
end

disp('Finished');
