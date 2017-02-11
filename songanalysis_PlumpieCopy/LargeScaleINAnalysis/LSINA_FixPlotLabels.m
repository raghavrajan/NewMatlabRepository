function [] = LSINA_FixPlotLabels(PlotHandle, PlotFontSize, varargin)

% A custom program to format all my plots with simple one line commands
% I decided to use parameter name, value pairs to do this so that I can be
% flexible about the things I put in to the program)


figure(PlotHandle);
set(gcf, 'Color', 'w');
set(gcf, 'Position', [100 100 900 450]);
set(gca, 'FontSize', PlotFontSize);

if (nargin > 2)
    NumParameters = (nargin/2) - 1;

    for i = 1:NumParameters,
        Parameter = varargin{(i-1)*2 + 1};
        Value = varargin{i*2};
        
        switch Parameter
            case 'Fig_Position'
                set(gcf, 'Position', Value);
            
            case 'Fig_Axis_Position'
                set(gca, 'Position', Value);
                
            case 'XAxis_Label'
                xlabel(Value, 'FontSize', PlotFontSize);
                
            case 'Fig_XTickLabels'
                set(gca, 'XTick', 1:1:length(Value), 'XTickLabel', Value);
                
            case 'Fig_XTickLabel_Rotation'
                set(gca, 'XTickLabelRotation', Value);
                
            case 'YAxis_Label'
                ylabel(Value, 'FontSize', PlotFontSize);
                
            case 'Fig_Title'
                title(Value, 'FontSize', PlotFontSize);
                
            case 'Fig_Legend'
                legend(Value);
                
            case 'Fig_Save'
                saveas(gcf, Value);
                
            case 'Fig_Axis'
                axis(Value);
        end
    end
end