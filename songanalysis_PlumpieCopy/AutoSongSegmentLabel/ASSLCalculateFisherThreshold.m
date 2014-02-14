function [Threshold] = ASSLCalculateFisherThreshold(Data)

MaxValue = -10000;
Threshold = min(Data);

PlotVals = [];
for i = linspace(min(Data), max(Data), 50),
    Group1 = Data(find(Data < i));
    Group2 = Data(find(Data >= i));
    
    Group1Indices = Data < i;
    
    FisherCrit1 = 2*abs((mean(Group1) - mean(Group2)))/sqrt((var(Group1) + var(Group2)));
    
    if (~isempty(Group1) && ~isempty(Group2))
        FisherCrit = sqrt(var(Data))/sqrt((var(Group1) + var(Group2)));
    else
        FisherCrit = 1;
    end
    
    if (FisherCrit > MaxValue)
        MaxValue = FisherCrit;
        Threshold = i;
    end
    PlotVals = [PlotVals; [i FisherCrit FisherCrit1]];
end
GroupMeansVars = [mean(Data(find(Data < Threshold))) mean(Data(find(Data >= Threshold))) var(Data(find(Data < Threshold))) var(Data(find(Data >= Threshold)))]; 
Threshold(2) = mean(Data(find(Data < Threshold))) + 2*std(Data(find(Data < Threshold)));

%disp('Calculated threshold');

% figure;
% plotyy(PlotVals(:,1), PlotVals(:,2), PlotVals(:,1), PlotVals(:,3));
% hold on;
% temp = axis;
% plot([Threshold Threshold], [temp(3) temp(4)], 'k--');
% 
