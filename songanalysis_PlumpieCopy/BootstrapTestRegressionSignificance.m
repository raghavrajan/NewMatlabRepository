function [FixedX_BootStrapCoeffs, Random_BootStrapCoeffs, RandomProb, Rsq] = BootstrapTestRegressionSignificance(Coeffs, XData, YData, FitType)

XData = XData(:);
YData = YData(:);

switch FitType
    case 'linear'
        YPred = Coeffs(2)*XData + Coeffs(1);
        

end

Dev = YData - mean(YData);          % deviations - measure of spread
SST = sum(Dev.^2);          % total variation to be accounted for
SSE = sum((YData(:) - YPred(:)).^2);

Rsq = 1 - SSE/SST;          % percent of error explained

YPred = YPred(:);

Error = YPred - YData;

BootStrapNumTrials = 2000;

rng('default');

% Fixed-x resampling
for i = 1:BootStrapNumTrials,
    BootStrapError = randsample(Error, length(Error), 'true');
    BootStrapYData = YPred + BootStrapError;
    
    switch FitType
        case 'linear'
            BootStrapCoeffs(i,:) = robustfit(XData, BootStrapYData);
    end
end

FixedX_BootStrapCoeffs = prctile(BootStrapCoeffs, [2.5 50 97.5]);

% Random-x resampling
clear BootStrapCoeffs
for i = 1:BootStrapNumTrials,
    BootStrapYData = randsample(YData, length(YData), 'true');
    switch FitType
        case 'linear'
            BootStrapCoeffs(i,:) = robustfit(XData, BootStrapYData);
    end
end
Random_BootStrapCoeffs = prctile(BootStrapCoeffs, [2.5 50 97.5]);

for i = 1:length(Coeffs),
    if (Coeffs(i) < Random_BootStrapCoeffs(2,i))
        RandomProb(i) = 1 - (length(find(BootStrapCoeffs(:,i) > Coeffs(i)))/BootStrapNumTrials);
    else
        RandomProb(i) = 1 - (length(find(BootStrapCoeffs(:,i) < Coeffs(i)))/BootStrapNumTrials);
    end
end
