function [Rsq, F] = CalculateGoodnessofLinearFit(coeff, x, y)

ypred = polyval(coeff,x);   % predictions
dev = y - mean(y);          % deviations - measure of spread
SST = sum(dev.^2);          % total variation to be accounted for
resid = y - ypred;          % residuals - measure of mismatch
SSE = sum(resid.^2);        % variation NOT accounted for
Rsq = 1 - SSE/SST;          % percent of error explained

SSR = SST - SSE;              % the "ANOVA identity"
dfr = 1;
dfe = length(x) - 1 - dfr;    % degrees of freedom for error
MSE = SSE/dfe;                % mean-square error of residuals
MSR = SSR/dfr;                % mean-square error for regression
F = MSR/MSE;                  % f-statistic for regression
