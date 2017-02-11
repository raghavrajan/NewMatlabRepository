function [] = PlotSAPFeatsFromMatFile(MatFile, Feature, MaxSyllLen)

for i = 1:length(MatFile),
    DataFiles{i} = load(MatFile{i});
    DotIndex = find(MatFile{i} == '.');
    ExptDate{i} = MatFile{i}((DotIndex(1)+1):(DotIndex(2)-1));
end


TempFields = fieldnames(DataFiles{i}.DirBout);
if (strfind(Feature, 'PC'))
    for i = 2:9,
        for Files = 1:length(DataFiles),
            TempFet = [];
            for j = 1:length(DataFiles{Files}.DirBout.FileName),
                TempFetVal = eval(['cell2mat(DataFiles{Files}.DirBout.',TempFields{i},'{j})']);
                TempFet = [TempFet; TempFetVal'];
            end
            DirFeats{Files}(:,i-1) = TempFet;

            TempFet = [];
            for j = 1:length(DataFiles{Files}.UnDirBout.FileName),
                TempFetVal = eval(['cell2mat(DataFiles{Files}.UnDirBout.',TempFields{i},'{j})']);
                TempFet = [TempFet; TempFetVal'];
            end
            UnDirFeats{Files}(:,i-1) = TempFet;
        end
    end

    FeatureVect = UnDirFeats{1};
    FeatureVect = zscore(FeatureVect);
    [Coeff, Score, Latent] = princomp(FeatureVect);

    Figures = findobj('Type', 'figure');

    if (length(Figures) > 0)
        FigureIndex = max(Figures) + 1;
    else
        FigureIndex = 1;
    end

    TempAxis = [];
    figure(FigureIndex);
    annotation('textbox',[0.4 0.9 0.2 0.05], 'FontSize',12,'FontWeight','bold','String','Directed song','LineStyle','none');
    figure(FigureIndex + 1);
    annotation('textbox',[0.4 0.9 0.2 0.05], 'FontSize',12,'FontWeight','bold','String','Undirected song','LineStyle','none');

    for i = 1:length(DirFeats),
        figure(FigureIndex);
        if (length(DataFiles) > 1)
            subplot(length(DataFiles)/2, 2, i);
        end
        DirFeatVect{i} = zscore(DirFeats{i});
        plot(DirFeatVect{i}*Coeff(:,1), DirFeatVect{i}*Coeff(:,2), 'k+', 'MarkerSize', 2);
        axis tight;
        TempAxis(end+1, :) = axis;

        figure(FigureIndex + 1);
        if (length(DataFiles) > 1)
            subplot(length(DataFiles)/2, 2, i);
        end
        UnDirFeatVect{i} = zscore(UnDirFeats{i});
        plot(UnDirFeatVect{i}*Coeff(:,1), UnDirFeatVect{i}*Coeff(:,2), 'k+', 'MarkerSize', 2);
        axis tight;
        TempAxis(end+1, :) = axis;
    end

    FinalAxis = [min(TempAxis(:,1)) max(TempAxis(:,2)) min(TempAxis(:,3)) max(TempAxis(:,4))];
    FinalAxis = TempAxis(2,:);

    for i = 1:length(DirFeats),
        figure(FigureIndex);
        if (length(DataFiles) > 1)
            subplot(length(DataFiles)/2, 2, i);
        end
        axis(FinalAxis);
        Pos = get(gca, 'Position');
        annotation('textbox',[(Pos(1) + Pos(3)/2) (Pos(2) + Pos(4)) 0.125 0.05], 'FontSize',12,'FontWeight','bold','String',ExptDate{i},'LineStyle','none');

        figure(FigureIndex + 1);
        if (length(DataFiles) > 1)
            subplot(length(DataFiles)/2, 2, i);
        end
        axis(FinalAxis);
        Pos = get(gca, 'Position');
        annotation('textbox',[(Pos(1) + Pos(3)/2) (Pos(2) + Pos(4)) 0.125 0.05], 'FontSize',12,'FontWeight','bold','String',ExptDate{i},'LineStyle','none');
    end

    RowEdges = linspace(FinalAxis(1), FinalAxis(2), 150);
    ColEdges = linspace(FinalAxis(3), FinalAxis(4), 150);

    Width = 1.5;
    GaussianLen = 2;
    Fs1 = 1;
    XGauss = 1:1:(1 + round(2 * GaussianLen * Width * (Fs1)));
    XGauss = XGauss - (length(XGauss) + 1)/2;
    GaussWin = (1/((Width * Fs1) * sqrt(2 * pi))) * exp(-(XGauss.*XGauss)/(2 * (Width * Fs1) * (Width * Fs1)));

    FigureIndex = FigureIndex + 2;
    for i = 1:length(DirFeats),
        figure(FigureIndex);
        if (length(DataFiles) > 1)
            subplot(length(DataFiles)/2, 2, i);
        end
        Output = hist2d(DirFeatVect{i}*Coeff(:,1), DirFeatVect{i}*Coeff(:,2), RowEdges, ColEdges);
        Output = Output/(sum(sum(Output)));
        Output = conv2(GaussWin, GaussWin, Output, 'same');
        contourf((RowEdges(1:end-1) + (RowEdges(2)-RowEdges(1))/2), (ColEdges(1:end-1) + (ColEdges(2) - ColEdges(1))/2), Output', 'LineStyle', 'none');
        colorbar;
        Pos = get(gca, 'Position');
        annotation('textbox',[(Pos(1) + Pos(3)/2) (Pos(2) + Pos(4)) 0.125 0.05], 'FontSize',12,'FontWeight','bold','String',ExptDate{i},'LineStyle','none');    

        figure(FigureIndex + 1);
        if (length(DataFiles) > 1)
            subplot(length(DataFiles)/2, 2, i);
        end
        Output = hist2d(UnDirFeatVect{i}*Coeff(:,1), UnDirFeatVect{i}*Coeff(:,2), RowEdges, ColEdges);
        Output = Output/(sum(sum(Output)));
        Output = conv2(GaussWin, GaussWin, Output, 'same');
        contourf((RowEdges(1:end-1) + (RowEdges(2)-RowEdges(1))/2), (ColEdges(1:end-1) + (ColEdges(2) - ColEdges(1))/2), Output', 'LineStyle', 'none');
        colorbar;
        Pos = get(gca, 'Position');
        annotation('textbox',[(Pos(1) + Pos(3)/2) (Pos(2) + Pos(4)) 0.125 0.05], 'FontSize',12,'FontWeight','bold','String',ExptDate{i},'LineStyle','none');
    end        
else
    for i = 2:8,
        if (length(strfind(TempFields{i}, Feature)) > 0)
            MatchIndex = i;
            break;
        end
    end
    
    ColIndex = 1;
    for i = [MatchIndex 9],
        for Files = 1:length(DataFiles),
            TempFet = [];
            for j = 1:length(DataFiles{Files}.DirBout.FileName),
                TempFetVal = eval(['cell2mat(DataFiles{Files}.DirBout.',TempFields{i},'{j})']);
                TempFet = [TempFet; TempFetVal'];
            end
            DirFeats{Files}(:,ColIndex) = TempFet;

            TempFet = [];
            for j = 1:length(DataFiles{Files}.UnDirBout.FileName),
                TempFetVal = eval(['cell2mat(DataFiles{Files}.UnDirBout.',TempFields{i},'{j})']);
                TempFet = [TempFet; TempFetVal'];
            end
            UnDirFeats{Files}(:,ColIndex) = TempFet;
        end
        ColIndex = ColIndex + 1;
    end

    Figures = findobj('Type', 'figure');

    if (length(Figures) > 0)
        FigureIndex = max(Figures) + 1;
    else
        FigureIndex = 1;
    end

    TempAxis = [];
    figure(FigureIndex);
    annotation('textbox',[0.4 0.9 0.2 0.05], 'FontSize',12,'FontWeight','bold','String',['Directed song - ', Feature],'LineStyle','none');
    figure(FigureIndex + 1);
    annotation('textbox',[0.4 0.9 0.2 0.05], 'FontSize',12,'FontWeight','bold','String',['Undirected song - ', Feature],'LineStyle','none');

    for i = 1:length(DirFeats),
        figure(FigureIndex);
        if (length(DataFiles) > 1)
            subplot(length(DataFiles)/2, 2, i);
        end
        plot(DirFeats{i}(:,2), DirFeats{i}(:,1), 'k+', 'MarkerSize', 2);
        axis tight;
        TempAxis(end+1, :) = axis;

        figure(FigureIndex + 1);
        if (length(DataFiles) > 1)
            subplot(length(DataFiles)/2, 2, i);
        end
        plot(UnDirFeats{i}(:,2), UnDirFeats{i}(:,1), 'k+', 'MarkerSize', 2);
        axis tight;
        TempAxis(end+1, :) = axis;
    end

    FinalAxis = [min(TempAxis(:,1)) max(TempAxis(:,2)) min(TempAxis(:,3)) max(TempAxis(:,4))];
    FinalAxis = TempAxis(2,:);
    FinalAxis(2) = MaxSyllLen;
    
     for i = 1:length(DirFeats),
        figure(FigureIndex);
        if (length(DataFiles) > 1)
            subplot(length(DataFiles)/2, 2, i);
        end
        axis(FinalAxis);
        Pos = get(gca, 'Position');
        annotation('textbox',[(Pos(1) + Pos(3)/2) (Pos(2) + Pos(4)) 0.125 0.05], 'FontSize',12,'FontWeight','bold','String',ExptDate{i},'LineStyle','none');

        figure(FigureIndex + 1);
        if (length(DataFiles) > 1)
            subplot(length(DataFiles)/2, 2, i);
        end
        axis(FinalAxis);
        Pos = get(gca, 'Position');
        annotation('textbox',[(Pos(1) + Pos(3)/2) (Pos(2) + Pos(4)) 0.125 0.05], 'FontSize',12,'FontWeight','bold','String',ExptDate{i},'LineStyle','none');
    end

    RowEdges = linspace(FinalAxis(1), FinalAxis(2), 150);
    ColEdges = linspace(FinalAxis(3), FinalAxis(4), 150);

    Width = 1.5;
    GaussianLen = 2;
    Fs1 = 1;
    XGauss = 1:1:(1 + round(2 * GaussianLen * Width * (Fs1)));
    XGauss = XGauss - (length(XGauss) + 1)/2;
    GaussWin = (1/((Width * Fs1) * sqrt(2 * pi))) * exp(-(XGauss.*XGauss)/(2 * (Width * Fs1) * (Width * Fs1)));

    FigureIndex = FigureIndex + 2;
    for i = 1:length(DirFeats),
        figure(FigureIndex);
        if (length(DataFiles) > 1)
            subplot(length(DataFiles)/2, 2, i);
        end
        Output = hist2d(DirFeats{i}(:,2), DirFeats{i}(:,1), RowEdges, ColEdges);
        Output = Output/(sum(sum(Output)));
        Output = conv2(GaussWin, GaussWin, Output, 'same');
        contourf((RowEdges(1:end-1) + (RowEdges(2)-RowEdges(1))/2), (ColEdges(1:end-1) + (ColEdges(2) - ColEdges(1))/2), Output', 'LineStyle', 'none');
        colorbar;
        Pos = get(gca, 'Position');
        annotation('textbox',[(Pos(1) + Pos(3)/2) (Pos(2) + Pos(4)) 0.125 0.05], 'FontSize',12,'FontWeight','bold','String',ExptDate{i},'LineStyle','none');    

        figure(FigureIndex + 1);
        if (length(DataFiles) > 1)
            subplot(length(DataFiles)/2, 2, i);
        end
        Output = hist2d(UnDirFeats{i}(:,2), UnDirFeats{i}(:,1), RowEdges, ColEdges);
        Output = Output/(sum(sum(Output)));
        Output = conv2(GaussWin, GaussWin, Output, 'same');
        contourf((RowEdges(1:end-1) + (RowEdges(2)-RowEdges(1))/2), (ColEdges(1:end-1) + (ColEdges(2) - ColEdges(1))/2), Output', 'LineStyle', 'none');
        colorbar;
        Pos = get(gca, 'Position');
        annotation('textbox',[(Pos(1) + Pos(3)/2) (Pos(2) + Pos(4)) 0.125 0.05], 'FontSize',12,'FontWeight','bold','String',ExptDate{i},'LineStyle','none');
    end
end
figure(FigureIndex);
annotation('textbox',[0.4 0.9 0.2 0.05], 'FontSize',12,'FontWeight','bold','String','Directed song','LineStyle','none');
figure(FigureIndex + 1);
annotation('textbox',[0.4 0.9 0.2 0.05], 'FontSize',12,'FontWeight','bold','String','Undirected song','LineStyle','none');
disp('Finished');