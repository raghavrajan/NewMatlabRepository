function [Thresholds, ErrorRate, NonLinearErrorRate, OptimalThreshold] = MA_EstimateFalsePositiveFalseNegativePercentage(CorrectMatches_Threshold, AllowedError)

% Now to estimate false positive percentage and false negative percentage
% for each threshold

Total = sum(CorrectMatches_Threshold(:,3));

TruePositive = sum(CorrectMatches_Threshold(:,end));
TrueNegative = Total - TruePositive;

for i = 1:(size(CorrectMatches_Threshold, 1)-1),
    Thresholds(i) = CorrectMatches_Threshold(i,2);
    
    NumAboveThreshold = sum(CorrectMatches_Threshold((i+1):end, 3));
    NumActualCorrect = sum(CorrectMatches_Threshold((i+1):end, 6));
    FalsePositive = NumAboveThreshold - NumActualCorrect;
    
    FalseNegative = sum(CorrectMatches_Threshold(1:i, 6));
    
    ErrorRate(i) = 1 - ((TruePositive + TrueNegative)/(TruePositive + FalseNegative + TrueNegative + FalsePositive));
    NonLinearErrorRate(i) = 1 - (sqrt((TruePositive/(TruePositive + FalseNegative)) * (TrueNegative/(TrueNegative + FalsePositive)))); 
end

% Optimal threshold is the threshold where sum of false positive percentage
% and false negative percentage is minimum

switch AllowedError
    case 0
        [MinVal, MinIndex] = min(ErrorRate);
        
    otherwise
        MinIndex = find(ErrorRate <= AllowedError, 1, 'first');
        if (isempty(MinIndex))
            [MinVal, MinIndex] = min(ErrorRate);
        end
end
OptimalThreshold = Thresholds(MinIndex);
        