function [Threshold] = ASSLCalculateFisherThreshold(Data)

MaxValue = -10000;
Threshold = min(Data);

PlotVals = [];
for i = linspace(min(Data), max(Data), 100),
    Group1 = Data < i;
    Group2 = Data >= i;
    
%    FisherCrit1 = ((mean(Group1) - mean(Group2))^2)/((var(Group1) + var(Group2)));
    
    if (~isempty(Group1) && ~isempty(Group2))
        FisherCrit = sqrt(CalculateVarForFisherThreshold(Data))/sqrt((CalculateVarForFisherThreshold(Data(Group1)) + CalculateVarForFisherThreshold(Data(Group2))));
    else
        FisherCrit = 1;
    end
    
    if (FisherCrit > MaxValue)
        MaxValue = FisherCrit;
        Threshold = i;
    end
    PlotVals = [PlotVals; [i FisherCrit]];
end
% GroupMeansVars = [mean(Data(Data < Threshold)) mean(Data(Data >= Threshold)) var(Data(Data < Threshold)) var(Data(Data >= Threshold))]; 
Threshold(2) = mean(Data(Data < Threshold)) + 2*std(Data(Data < Threshold));
if (Threshold(2) > Threshold(1))
    Threshold(2) = Threshold(1) - 1;
end
%disp('Calculated threshold');

% figure;
% plotyy(PlotVals(:,1), PlotVals(:,2), PlotVals(:,1), PlotVals(:,3));
% hold on;
% temp = axis;
% plot([Threshold Threshold], [temp(3) temp(4)], 'k--');
% 
