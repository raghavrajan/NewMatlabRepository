function [] = PolarPlotsForGayathri(InputFile, FilledOrNot, varargin)

Fid = fopen(InputFile, 'r');
Text = textscan(Fid, '%s', 'DeLimiter', '\n');
Text = Text{1};
fclose(Fid);

if (nargin > 2)
    StarPointResidueNames = varargin{1};
    StarPointAngles = varargin{2};
end

for i = 1:length(Text),
    LineText = textscan(Text{i}, '%s', 'DeLimiter', ',');
    LineText = LineText{1};
    Theta(i) = str2double(LineText{4});
    Phi(i) = str2double(LineText{6});
    AtomName{i} = LineText{1};
    ResidueName{i} = LineText{2};
end

Phi = 2*pi * ((Phi + 180)/360);

SizeofMarker = 4;
LineThickness = 2;
Os = strmatch('O', AtomName, 'exact');
Legend = {};
if (strcmp('filled', FilledOrNot) == 1)
    PolarPlotHandles(1) = polarplot(Phi(Os), Theta(Os), 'bo', 'MarkerSize', SizeofMarker, 'MarkerFaceColor', 'b');
else
    PolarPlotHandles(1) = polarplot(Phi(Os), Theta(Os), 'bo', 'MarkerSize', SizeofMarker);
end
hold on;
% if (length(Os) < 5)
%     polarplot([0 mean(Phi(Os))], [0 mean(Theta(Os))], 'b:', 'LineWidth', LineThickness-0.3);
% else
%     Phi_X = Theta(Os).*cos(Phi(Os));
%     Phi_Y = Theta(Os).*sin(Phi(Os));
%     MeanPhi_X = mean(Phi_X);
%     MeanPhi_Y = mean(Phi_Y);
%     MeanTheta = sqrt(MeanPhi_X.^2 + MeanPhi_Y.^2);
%     MeanPhi = atan(MeanPhi_Y/MeanPhi_X);
%     polarplot([0 MeanPhi], [0 MeanTheta], 'b', 'LineWidth', LineThickness);
% end
        
Legend{end+1} = 'Os';

NotOs = setdiff(1:1:length(AtomName), Os);

SpecificResidues = {'GLU' 'ASP' 'SER' 'THR' 'ASN' 'GLN' 'TYR'};

ColourForResidues = [0 0 0; 0.8 0.8 0.8; 0 1 0; 1 0 1; 1 0.87 0.37; 1 0.5 0; 1 0 0.5];

UniqueResidues = unique(ResidueName);
for i = 1:length(UniqueResidues),
    ResidueIndices = strmatch(UniqueResidues{i}, ResidueName(NotOs), 'exact');
    ResidueIndices = NotOs(ResidueIndices);
    if (~isempty(ResidueIndices))
        ColourIndex = strmatch(UniqueResidues{i}, SpecificResidues, 'exact');
        if (strcmp('filled', FilledOrNot) == 1)
            PolarPlotHandles(end+1) = polarplot(Phi(ResidueIndices), Theta(ResidueIndices), 'ko', 'Color', ColourForResidues(ColourIndex,:), 'MarkerSize', SizeofMarker, 'MarkerFaceColor', ColourForResidues(ColourIndex,:));
        else
            PolarPlotHandles(end+1) = polarplot(Phi(ResidueIndices), Theta(ResidueIndices), 'ko', 'Color', ColourForResidues(ColourIndex,:), 'MarkerSize', SizeofMarker);
        end
        
        Phi_X = Theta(ResidueIndices).*cos(Phi(ResidueIndices));
        Phi_Y = Theta(ResidueIndices).*sin(Phi(ResidueIndices));
        MeanPhi_X = mean(Phi_X);
        MeanPhi_Y = mean(Phi_Y);
        MeanTheta = sqrt(MeanPhi_X.^2 + MeanPhi_Y.^2);
        MeanPhi = atan(MeanPhi_Y/MeanPhi_X);
    
%         if (length(ResidueIndices) < 5)
%             polarplot([0 MeanPhi], [0 MeanTheta], 'k:', 'Color', ColourForResidues(ColourIndex,:), 'LineWidth', LineThickness-0.3);
%         else
%             polarplot([0 MeanPhi], [0 MeanTheta], 'k', 'Color', ColourForResidues(ColourIndex,:), 'LineWidth', LineThickness-0.3);
%         end
        Legend{end+1} = UniqueResidues{i};
    end
end

if (nargin > 2)
    StarPointAngles(:,1) = 2*pi * ((StarPointAngles(:,1) + 180)/360);
    for i = 1:length(StarPointResidueNames),
        ResidueIndices = strmatch(StarPointResidueNames{i}, SpecificResidues, 'exact');
        if (~isempty(ResidueIndices))
            polarplot(StarPointAngles(i,1), StarPointAngles(i,2), 'k*', 'Color', ColourForResidues(ResidueIndices,:), 'MarkerSize', SizeofMarker*2, 'MarkerFaceColor', ColourForResidues(ResidueIndices,:));
        else
            polarplot(StarPointAngles(i,1), StarPointAngles(i,2), 'k*', 'Color', 'b', 'MarkerSize', SizeofMarker*2, 'MarkerFaceColor', 'b');
        end
    end
end


legend(PolarPlotHandles, Legend);
set(gca, 'FontSize', 15, 'FontWeight', 'bold');
set(gcf, 'PaperPositionMode', 'auto');
set(gcf, 'Position', [680 403 750 600]);
print(gcf, [InputFile, '.jpg'], '-djpeg', '-r600');
print(gcf, [InputFile, '.eps'], '-depsc2', '-r600');
disp('FInished');