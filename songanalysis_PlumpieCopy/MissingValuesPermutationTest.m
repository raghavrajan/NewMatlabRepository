function [RobustCoeffs, RandomPosPValue, SamePosPValue] = MissingValuesPermutationTest(XValues, Data)

% First generate the different missing value randomized datasets keeping
% the values in a row always within that row itself.
% Two ways of doing this:
% 1) Keep missing value positions always the same
% 2) Randomise missing value positions

for i = 1:10000,
    MissingPosSameDataSet = ones(size(Data))*NaN;
    MissingPosRandomDataSet = ones(size(Data))*NaN;
    % Type 1 with missing value positions the same
    for j = 1:size(Data, 1),
        NonNaNValueIndices = find(~isnan(Data(j,:)));
        
        NonNaNValues = Data(j,NonNaNValueIndices);
        
        MissingPosSameDataSet(j,NonNaNValueIndices) = randsample(NonNaNValues, length(NonNaNValues), 'false');
        
        MissingPosRandomDataSet(j,randsample(1:1:size(Data,2), length(NonNaNValues), 'false')) = randsample(NonNaNValues, length(NonNaNValues), 'false');
    end
    MissingPosSameCoeffs(i,:) = robustfit(XValues(:), MissingPosSameDataSet(:));
    MissingPosRandomCoeffs(i,:) = robustfit(XValues(:), MissingPosRandomDataSet(:));
end

RobustCoeffs = robustfit(XValues(:), Data(:));

% figure;
% plot(XValues', Data', 'ko-', 'Color', [0.75 0.75 0.75]);
% hold on;
% for i = 1:size(Data,2),
%     errorbar(mean(XValues(:,i)), nanmean(Data(:,i)), nanstd(Data(:,i))/sqrt(length(find(~isnan(Data(:,i))))), 'ks', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
% end
% plot(mean(XValues), polyval(flipud(RobustCoeffs), mean(XValues)), 'k', 'LineWidth', 1.5);

% Plot mean and 95% confidence intervals for the two random sets in blue
% and red
% plot(mean(XValues), polyval(fliplr(mean(MissingPosSameCoeffs)), mean(XValues)), 'b', 'LineWidth', 1);
% plot(mean(XValues), polyval(fliplr(prctile(MissingPosSameCoeffs, 2.5)), mean(XValues)), 'b--', 'LineWidth', 1);
% plot(mean(XValues), polyval(fliplr(prctile(MissingPosSameCoeffs, 97.5)), mean(XValues)), 'b--', 'LineWidth', 1);
% 
% plot(mean(XValues), polyval(fliplr(mean(MissingPosRandomCoeffs)), mean(XValues)), 'r', 'LineWidth', 1);
% plot(mean(XValues), polyval(fliplr(prctile(MissingPosRandomCoeffs, 2.5)), mean(XValues)), 'r--', 'LineWidth', 1);
% plot(mean(XValues), polyval(fliplr(prctile(MissingPosRandomCoeffs, 97.5)), mean(XValues)), 'r--', 'LineWidth', 1);

disp(RobustCoeffs);
disp(prctile(MissingPosRandomCoeffs, [2.5 97.5]));
disp(prctile(MissingPosSameCoeffs, [2.5 97.5]));

for i = 1:length(RobustCoeffs),
    if (RobustCoeffs(i) >= prctile(MissingPosRandomCoeffs(:,i), 50))
        disp(['Random position p-value = ', num2str(1 - length(find(MissingPosRandomCoeffs(:,i) <= RobustCoeffs(i)))/length(MissingPosRandomCoeffs))]);
        RandomPosPValue(i) = 1 - (length(find(MissingPosRandomCoeffs(:,i) <= RobustCoeffs(i)))/length(MissingPosRandomCoeffs));
    else
        disp(['Random position p-value = ', num2str(1 - length(find(MissingPosRandomCoeffs(:,i) >= RobustCoeffs(i)))/length(MissingPosRandomCoeffs))]);
        RandomPosPValue(i) = 1 - (length(find(MissingPosRandomCoeffs(:,i) >= RobustCoeffs(i)))/length(MissingPosRandomCoeffs));
    end

    if (RobustCoeffs(i) >= prctile(MissingPosSameCoeffs(:,i), 50))
        disp(['Same position p-value = ', num2str(1 - length(find(MissingPosSameCoeffs(:,i) <= RobustCoeffs(i)))/length(MissingPosSameCoeffs))]);
        SamePosPValue(i) = 1 - length(find(MissingPosSameCoeffs(:,i) <= RobustCoeffs(i)))/length(MissingPosSameCoeffs);
    else
        disp(['Same position p-value = ', num2str(1 - length(find(MissingPosSameCoeffs(:,i) >= RobustCoeffs(i)))/length(MissingPosSameCoeffs))]);
        SamePosPValue(i) = 1 - length(find(MissingPosSameCoeffs(:,i) >= RobustCoeffs(i)))/length(MissingPosSameCoeffs);
    end
end
disp('Finshed');