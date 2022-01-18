function [] = SimulateIQR_NumRelationship(NumIndependents, MeanDifference, STDDifference, FractionDifferent)

% This is to check whether iqr varies with number of samples and how it
% varies.

%% Variables
for i = 1:NumIndependents,
    NumSessions(i) = ceil((rand * 5) + 3);
end

%% Simulating the session data
% No difference in mean between sessions
SessionMean = 4;
SessionSigma = 0.5;

for i = 1:NumIndependents,
    NumSessionsDifferent = round(FractionDifferent * NumSessions(i));
    for j = 1:NumSessions(i),
        SessionNumBouts{i}(j) = ceil(rand*10) + 3;
        if (j <= NumSessionsDifferent)
            SessionValues{i}{j} = normrnd(SessionMean + normrnd(MeanDifference, STDDifference), SessionSigma, SessionNumBouts{i}(j), 1);
        else
            SessionValues{i}{j} = normrnd(SessionMean, SessionSigma, SessionNumBouts{i}(j), 1);
        end
    end
end

%% Now measure median and iqr for sessions and global
NumResamples = 10000;
for i = 1:NumIndependents,
    AllSessionVals = cell2mat(SessionValues{i}(:));
    for j = 1:NumSessions(i),
        SessionMedian{i}(j) = median(SessionValues{i}{j});
        SessionIQR{i}(j) = iqr(SessionValues{i}{j});
        
        TempMedian = ones(NumResamples, 1)*NaN;
        TempIQR = ones(NumResamples, 1)*NaN;
        
        for k = 1:10000,
            TempIndices = randperm(length(cell2mat(SessionValues{i}(:))));
            TempIndices = TempIndices(1:length(SessionValues{i}{j}));
            TempMedian(k) = median(AllSessionVals(TempIndices));
            TempIQR(k) = iqr(AllSessionVals(TempIndices));
        end
        SessionMedian_CIs{i}{j} = prctile(TempMedian, [2.5 97.5]);
        SessionIQR_CIs{i}{j} = prctile(TempIQR, [2.5 97.5]);
    end
    GlobalMedian(i) = median(AllSessionVals);
    GlobalIQR(i) = iqr(AllSessionVals);
end

%% Now to plot data
figure;
p = panel();
p.pack(NumIndependents, 1);
% First plot data
for i = 1:NumIndependents,
    p(i,1).select();
    hold on;
    for j = 1:NumSessions(i),
        ColourBoxPlot(j, SessionValues{i}{j}, 'k', 0.4, 'filled');
        plot(j+[-0.4 0.4], ones(1,2)*SessionMedian_CIs{i}{j}(1), 'r');
        plot(j+[-0.4 0.4], ones(1,2)*SessionMedian_CIs{i}{j}(2), 'r');
    end
    if (i == 1)
        title({['Fraction of sessions with different mean = ', num2str(FractionDifferent)]; ['Session mean difference = ', num2str(MeanDifference)]; ['Session sigma difference = ', num2str(STDDifference)]});
    end
end
xlabel('Session Index');
ylabel('Measurement');
p.fontsize = 16;
set(gcf, 'Position', [150 94 1000 900]);
p.margintop = 20;

figure;
p = panel();
p.pack(2,2);
% For each independent, I should plot the individual sessions that are
% outside the 95% CI as filled circles and the rest as open circles
for i = 1:NumIndependents,
    for j = 1:NumSessions(i),
        p(1,1).select();
        hold on;
        if ((SessionMedian{i}(j) >= SessionMedian_CIs{i}{j}(1)) & (SessionMedian{i}(j) <= SessionMedian_CIs{i}{j}(2)))
            plot(i, SessionMedian{i}(j), 'ko');
        else
            plot(i, SessionMedian{i}(j), 'ko', 'MarkerFaceColor', 'k');
        end
        p(1,2).select();
        hold on;
        if ((SessionIQR{i}(j) >= SessionIQR_CIs{i}{j}(1)) & (SessionIQR{i}(j) <= SessionIQR_CIs{i}{j}(2)))
            plot(i, SessionIQR{i}(j), 'ko');
        else
            plot(i, SessionIQR{i}(j), 'ko', 'MarkerFaceColor', 'k');
        end
    end
end
p(1,1).select();
axis tight;
Temp = axis;
Temp = [0.5 NumIndependents+0.5 0.98*Temp(3) 1.02*Temp(4)];
axis(Temp);
ylabel('Session medians');
xlabel('Independents');
title({['Fraction of sessions with different mean = ', num2str(FractionDifferent)]; ['Session mean difference = ', num2str(MeanDifference)]; ['Session sigma difference = ', num2str(STDDifference)]});

p(1,2).select();
axis tight;
Temp = axis;
Temp = [0.5 NumIndependents+0.5 0.98*Temp(3) 1.02*Temp(4)];
axis(Temp);
ylabel('Session IQRs');
xlabel('Independents');
title({['Fraction of sessions with different mean = ', num2str(FractionDifferent)]; ['Session mean difference = ', num2str(MeanDifference)]; ['Session sigma difference = ', num2str(STDDifference)]});

% Now plot means for medians and iqrs as paired data
p(2,1).select();
hold on;
errorbar(0.8, mean(cellfun(@mean, SessionMedian)), std(cellfun(@mean, SessionMedian))/sqrt(length(SessionMedian)), 'ks', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
plot(repmat([1 2], length(SessionMedian), 1)', [cellfun(@mean, SessionMedian); GlobalMedian(:)'], 'ko-', 'Color', [0.5 0.5 0.5]);
errorbar(2.2, mean(GlobalMedian), std(GlobalMedian)/sqrt(length(GlobalMedian)), 'ks', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
axis tight;
Temp = axis;
Temp = [0.3 2.7 0.98*Temp(3) 1.02*Temp(4)];
axis(Temp);
set(gca, 'XTick', [1 2], 'XTickLabel', {'Session means' 'Global means'});
ylabel('Median');

p(2,2).select();
hold on;
errorbar(0.8, mean(cellfun(@mean, SessionIQR)), std(cellfun(@mean, SessionIQR))/sqrt(length(SessionIQR)), 'ks', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
plot(repmat([1 2], length(SessionIQR), 1)', [cellfun(@mean, SessionIQR); GlobalIQR(:)'], 'ko-', 'Color', [0.5 0.5 0.5]);
errorbar(2.2, mean(GlobalIQR), std(GlobalIQR)/sqrt(length(GlobalIQR)), 'ks', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
axis tight;
Temp = axis;
Temp = [0.3 2.7 0.98*Temp(3) 1.02*Temp(4)];
axis(Temp);
set(gca, 'XTick', [1 2], 'XTickLabel', {'Session means' 'Global means'});
ylabel('IQRs');
p.margintop = 20;
set(gcf, 'Position', [680 195 1050 900]);
p.fontsize = 14;

disp('Finished');