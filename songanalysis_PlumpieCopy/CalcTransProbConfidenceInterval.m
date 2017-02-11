function [CI] = CalcTransProbConfidenceInterval(NumberofSamples, TransProb, NumRepetitions, Alpha)

Label = [repmat('bb', round(TransProb * NumberofSamples), 1); repmat('bc', round((1-TransProb) * NumberofSamples), 1)];

for i = 1:NumRepetitions,
    Indices = ceil(rand(NumberofSamples, 1)*NumberofSamples);
    TrialLabels = Label(Indices,:);
    TrialLabels = mat2cell(TrialLabels, ones(NumberofSamples, 1), 2);
    TransProb(i) = length(find(cellfun(@length, strfind(TrialLabels, 'bb'))));
end

TransProb = sort(TransProb);

CI = [TransProb(round(Alpha/2*NumRepetitions)) TransProb(round((1-Alpha/2)*NumRepetitions))]/NumberofSamples;

figure;
plot([0:1:NumberofSamples]/NumberofSamples, histc(TransProb, [0:1:NumberofSamples])/sum(histc(TransProb, [0:1:NumberofSamples])), 'k');
xlabel('Transition Probability', 'FontSize', 16);
ylabel('Fraction of trials', 'FontSize', 16);
