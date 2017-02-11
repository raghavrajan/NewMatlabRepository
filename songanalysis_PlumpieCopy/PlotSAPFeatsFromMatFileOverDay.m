function [] = PlotSAPFeatsFromMatFileOverDay(MatFile, Feature, MaxSyllLen, Date)

TimeInts = [8.5 10.5; 10.5 12.5; 12.5 14.5; 14.5 16.5; 16.5 18.5];

for i = 1:length(MatFile),
    DataFiles{i} = load(MatFile{i});
    DotIndex = find(MatFile{i} == '.');
    BirdName = MatFile{i}(1:DotIndex(1)-1);
    ExptDate{i} = MatFile{i}((DotIndex(1)+1):(DotIndex(2)-1));
    if (strfind(ExptDate{i}, Date))
        DateIndex = i;
    end
end

FileName = [];
FileTime = [];
TimeFileFid = fopen([MatFile{DateIndex}(1:DotIndex(2)), 'times'], 'r');
Tline = fgetl(TimeFileFid);
FileInd = 1;
while (length(strfind(Tline, 'Total')) < 1)
    if (length(Tline) > 0)
        if ((strfind(Tline, BirdName)) & (strfind(Tline, 'Recorded')))
            ColonIndex = find(Tline == ':');
            CommaIndex = find(Tline == ',');
            FileName{FileInd} = Tline(ColonIndex(1) + 2:CommaIndex(1) - 1);
            TempFileTime = str2double(Tline(ColonIndex(3) - 2:ColonIndex(3) - 1)) + str2double(Tline(ColonIndex(3)+1:ColonIndex(4)-1))/60 + str2double(Tline(ColonIndex(4)+1:ColonIndex(4)+2))/3600;
            FileTime = [FileTime; TempFileTime];
            FileInd = FileInd + 1;
        end
    end
    Tline = fgetl(TimeFileFid);
end
fclose(TimeFileFid);

TempFields = fieldnames(DataFiles{i}.DirBout);
if (strfind(Feature, 'PC'))
    for i = 2:9,
        for Files = DateIndex,
            TempFet = [];
            for j = 1:length(DataFiles{Files}.DirBout.FileName),
                TempFetVal = eval(['cell2mat(DataFiles{Files}.DirBout.',TempFields{i},'{j})']);
                TempFet = [TempFet; TempFetVal'];
            end
            DirFeats{1}(:,i-1) = TempFet;

            TempFet = [];
            for j = 1:length(DataFiles{Files}.UnDirBout.FileName),
                TempFetVal = eval(['cell2mat(DataFiles{Files}.UnDirBout.',TempFields{i},'{j})']);
                TempFet = [TempFet; TempFetVal'];
            end
            UnDirFeats{1}(:,i-1) = TempFet;
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
        for Files = DateIndex,
            TempFet = [];
            TempFileTime = [];
            for j = 1:length(DataFiles{Files}.DirBout.FileName),
                TempFetVal = eval(['cell2mat(DataFiles{Files}.DirBout.',TempFields{i},'{j})']);
                TempFet = [TempFet; TempFetVal'];
                TempMatches = strfind(FileName, DataFiles{Files}.DirBout.FileName{j}{1});
                FileTimeIndex = find(cellfun(@length, TempMatches));
                TempFileTime = [TempFileTime; ones(size(TempFetVal'))*FileTime(FileTimeIndex)];
            end
            DirFeats{1}(:,ColIndex) = TempFet;
            DirFeatTimes = TempFileTime;
            
            TempFet = [];
            TempFileTime = [];
            for j = 1:length(DataFiles{Files}.UnDirBout.FileName),
                TempFetVal = eval(['cell2mat(DataFiles{Files}.UnDirBout.',TempFields{i},'{j})']);
                TempFet = [TempFet; TempFetVal'];
                TempMatches = strfind(FileName, DataFiles{Files}.UnDirBout.FileName{j}{1});
                FileTimeIndex = find(cellfun(@length, TempMatches));
                TempFileTime = [TempFileTime; ones(size(TempFetVal'))*FileTime(FileTimeIndex)];
            end
            UnDirFeats{1}(:,ColIndex) = TempFet;
            UnDirFeatTimes = TempFileTime;
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
    annotation('textbox',[0.4 0.9 0.2 0.05], 'FontSize',12,'FontWeight','bold','String',[BirdName, ' - ', Date, ' - ', Feature],'LineStyle','none');
    
    for i = 1:size(TimeInts, 1),
        figure(FigureIndex);
        subplot(size(TimeInts,1), 2, (i-1)*2 + 1);
        Indices = find((DirFeatTimes > TimeInts(i,1)) & (DirFeatTimes < TimeInts(i,2)));
        if (length(Indices) > 0)
            plot(DirFeats{1}(Indices,2), DirFeats{1}(Indices,1), 'k+', 'MarkerSize', 2);
            axis tight;
            TempAxis(end+1, :) = axis;
        end
        subplot(size(TimeInts,1), 2, i*2);
        Indices = find((UnDirFeatTimes > TimeInts(i,1)) & (UnDirFeatTimes < TimeInts(i,2)));
        if (length(Indices) > 0)
            plot(UnDirFeats{1}(Indices,2), UnDirFeats{1}(Indices,1), 'k+', 'MarkerSize', 2);
            axis tight;
            TempAxis(end+1, :) = axis;
        end
    end

    FinalAxis = [min(TempAxis(:,1)) max(TempAxis(:,2)) min(TempAxis(:,3)) max(TempAxis(:,4))];
    FinalAxis = TempAxis(2,:);
    FinalAxis(2) = MaxSyllLen;
    
    for i = 1:size(TimeInts,1),
        figure(FigureIndex);
        subplot(size(TimeInts, 1), 2, (i-1)*2 + 1);
        axis(FinalAxis);
        Pos = get(gca, 'Position');
        annotation('textbox',[(Pos(1) + Pos(3)/2) (Pos(2) + Pos(4)) 0.125 0.05], 'FontSize',12,'FontWeight','bold','String',[num2str(TimeInts(i,1)), '-', num2str(TimeInts(i,2))],'LineStyle','none');

        subplot(size(TimeInts, 1), 2, i*2);
        axis(FinalAxis);
        Pos = get(gca, 'Position');
        annotation('textbox',[(Pos(1) + Pos(3)/2) (Pos(2) + Pos(4)) 0.125 0.05], 'FontSize',12,'FontWeight','bold','String',[num2str(TimeInts(i,1)), '-', num2str(TimeInts(i,2))],'LineStyle','none');
    end

    RowEdges = linspace(FinalAxis(1), FinalAxis(2), 150);
    ColEdges = linspace(FinalAxis(3), FinalAxis(4), 150);

    Width = 1.5;
    GaussianLen = 2;
    Fs1 = 1;
    XGauss = 1:1:(1 + round(2 * GaussianLen * Width * (Fs1)));
    XGauss = XGauss - (length(XGauss) + 1)/2;
    GaussWin = (1/((Width * Fs1) * sqrt(2 * pi))) * exp(-(XGauss.*XGauss)/(2 * (Width * Fs1) * (Width * Fs1)));

    FigureIndex = FigureIndex + 1;
    for i = 1:size(TimeInts,1),
        figure(FigureIndex);
        subplot(size(TimeInts,1), 2, (i-1)*2 + 1);
        Indices = find((DirFeatTimes > TimeInts(i,1)) & (DirFeatTimes < TimeInts(i,2)));
        Output = hist2d(DirFeats{1}(Indices,2), DirFeats{1}(Indices,1), RowEdges, ColEdges);
        Output = Output/(sum(sum(Output)));
        Output = conv2(GaussWin, GaussWin, Output, 'same');
        contourf((RowEdges(1:end-1) + (RowEdges(2)-RowEdges(1))/2), (ColEdges(1:end-1) + (ColEdges(2) - ColEdges(1))/2), Output', 'LineStyle', 'none');
        colorbar;
        Pos = get(gca, 'Position');
        annotation('textbox',[(Pos(1) + Pos(3)/2) (Pos(2) + Pos(4)) 0.125 0.05], 'FontSize',12,'FontWeight','bold','String',[num2str(TimeInts(i,1)), '-', num2str(TimeInts(i,2))],'LineStyle','none');    

        subplot(size(TimeInts,1), 2, i*2);
        Indices = find((UnDirFeatTimes > TimeInts(i,1)) & (UnDirFeatTimes < TimeInts(i,2)));
        Output = hist2d(UnDirFeats{1}(Indices,2), UnDirFeats{1}(Indices,1), RowEdges, ColEdges);
        Output = Output/(sum(sum(Output)));
        Output = conv2(GaussWin, GaussWin, Output, 'same');
        contourf((RowEdges(1:end-1) + (RowEdges(2)-RowEdges(1))/2), (ColEdges(1:end-1) + (ColEdges(2) - ColEdges(1))/2), Output', 'LineStyle', 'none');
        colorbar;
        Pos = get(gca, 'Position');
        annotation('textbox',[(Pos(1) + Pos(3)/2) (Pos(2) + Pos(4)) 0.125 0.05], 'FontSize',12,'FontWeight','bold','String',[num2str(TimeInts(i,1)), '-', num2str(TimeInts(i,2))],'LineStyle','none');
    end
end
figure(FigureIndex);
annotation('textbox',[0.4 0.9 0.2 0.05], 'FontSize',12,'FontWeight','bold','String',[BirdName, ' - ', Date, ' - ', Feature], 'LineStyle','none');
disp('Finished');
