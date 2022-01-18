function [] = MarkerSize_AxisLimits_Relation()

% This is to a function to understand the relation between marker size
% (circle) and the spacing between two markers that will make them distinct
% - as a function of the limits of the axis.

MarkerSymbol = '.';
MarkerSizes = [2 4 6 8];
InterMarkerDistances = 0.01:0.01:0.2;

figure;
set(gcf, 'Color', 'w');
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [1 16 9 21]);
p = panel();
p.pack({1/4 1/4 1/4 1/4});


% First for a set of 3 groups (the axis limits will be 0.5 to 3.5)
% Vary distance between two points by 0.02 to 0.2;
for i = 1:4,
    p(i).select();
    hold on;
    YTickLabelString = [];
    for j = 1:length(InterMarkerDistances),
        plot(1,j,['k', MarkerSymbol], 'MarkerSize', MarkerSizes(i), 'MarkerFaceColor', 'k');
        plot(1+InterMarkerDistances(j), j, ['k', MarkerSymbol], 'MarkerSize', MarkerSizes(i), 'MarkerFaceColor', 'k');
        YTickLabelString{j} = InterMarkerDistances(j);
    end
    axis([0.5 5.5 0.5 j+0.5]);
    ylabel('Inter marker distances');
    set(gca, 'YTick', 1:3:j, 'YTickLabel', YTickLabelString(1:3:end));
    title(['Marker size = ', num2str(MarkerSizes(i))]);
end

p.fontsize = 8;
p.marginleft = 10;
