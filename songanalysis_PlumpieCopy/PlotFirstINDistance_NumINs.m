function [] = PlotFirstINDistance_NumINs(IntroNoteFeatureAnalysis_Result)

figure;
set(gcf, 'Color', 'w');
hold on;
FirstINsDistanceToLast = [];
FirstINsDistanceToLast = [FirstINsDistanceToLast; IntroNoteFeatureAnalysis_Result.FirstINsDistanceToLast];
plot(IntroNoteFeatureAnalysis_Result.FirstINsDistanceToLast(:,1), IntroNoteFeatureAnalysis_Result.FirstINsDistanceToLast(:,2), 'ko', 'MarkerSize', 4);

for i = min(FirstINsDistanceToLast(:,1)):1:max(FirstINsDistanceToLast(:,1)),
    errorbar(i+0.2, mean(FirstINsDistanceToLast(find(FirstINsDistanceToLast(:,1) == i),2)), std(FirstINsDistanceToLast(find(FirstINsDistanceToLast(:,1) == i),2)), 'ks-', 'MarkerFaceColor', 'k');
end

[r, p] = corrcoef(FirstINsDistanceToLast(:,1), FirstINsDistanceToLast(:,2));
disp(['Correlation between first IN distance to last and number of INs is r = ', num2str(r(1,2)), ' and p = ', num2str(p(1,2))]);

axis tight;
temp = axis;
temp = [(temp(1)-0.2) (temp(2)+0.2) 0 22];
axis(temp);
set(gca, 'YTick', [0:5:20]);
set(gca, 'XTick', [min(FirstINsDistanceToLast(:,1)):1:max(FirstINsDistanceToLast(:,1))]);
set(gca, 'FontSize', 16, 'FontName', 'Times New Roman');
xlabel('# of INs', 'FontSize', 16);
ylabel('Distance of First IN from last IN position', 'FontSize', 16);
