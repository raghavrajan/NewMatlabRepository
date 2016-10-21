function [] = MA_ProcessSyllableAnalysisData(Parameters, Results)

Colors = 'rgbcmy';
Symbols = 'o+*<>';

for i = 1:length(Parameters.SyllableTemplate.SyllableTemplates),
    SyllLabels(i) = Parameters.SyllableTemplate.SyllableTemplates{i}{1}.MotifTemplate(1).Label;
end

for i = 1:Parameters.NoPreDays,
    for j = 1:length(Parameters.SyllablePreUnDirResults{i}{1}),
        PreUnDirMatchTimes{i}{j} = [];
        PreUnDirMatchLabels{i}{j} = [];
        PreUnDirPlotLabels{i}{j} = [];
        for k = 1:length(SyllLabels),
            if (~isempty(Parameters.SyllablePreUnDirResults{i}{k}{j}))
                Indices = find(Parameters.SyllablePreUnDirResults{i}{k}{j}(:,1) >= Results.ActualSyllableTemplateMatchThreshold{k});
                PreUnDirMatchTimes{i}{j} = [PreUnDirMatchTimes{i}{j}; Parameters.SyllablePreUnDirResults{i}{k}{j}(Indices,2)];
                PreUnDirMatchLabels{i}{j} = [PreUnDirMatchLabels{i}{j}; char(ones(size(Indices))*SyllLabels(k))];
                PreUnDirPlotLabels{i}{j} = [PreUnDirPlotLabels{i}{j}; (ones(size(Indices))*k)];
            end
        end
        if (~isempty(PreUnDirMatchLabels{i}{j}))
            [SortedVals, SortedIndices] = sort(PreUnDirMatchTimes{i}{j});
            PreUnDirMatchTimes{i}{j} = PreUnDirMatchTimes{i}{j}(SortedIndices);
            PreUnDirMatchLabels{i}{j} = PreUnDirMatchLabels{i}{j}(SortedIndices);
            PreUnDirPlotLabels{i}{j} = PreUnDirPlotLabels{i}{j}(SortedIndices);
        end
    end
end

for i = 1:Parameters.NoPostDays,
    for j = 1:length(Parameters.SyllablePostUnDirResults{i}{1}),
        PostUnDirMatchTimes{i}{j} = [];
        PostUnDirMatchLabels{i}{j} = [];
        PostUnDirPlotLabels{i}{j} = [];
        for k = 1:length(SyllLabels),
            if (~isempty(Parameters.SyllablePostUnDirResults{i}{k}{j}))
                Indices = find(Parameters.SyllablePostUnDirResults{i}{k}{j}(:,1) >= Results.ActualSyllableTemplateMatchThreshold{k});
                PostUnDirMatchTimes{i}{j} = [PostUnDirMatchTimes{i}{j}; Parameters.SyllablePostUnDirResults{i}{k}{j}(Indices,2)];
                PostUnDirMatchLabels{i}{j} = [PostUnDirMatchLabels{i}{j}; char(ones(size(Indices))*SyllLabels(k))];
                PostUnDirPlotLabels{i}{j} = [PostUnDirPlotLabels{i}{j}; (ones(size(Indices))*k)];
            end
        end
        if (~isempty(PostUnDirMatchLabels{i}{j}))
            [SortedVals, SortedIndices] = sort(PostUnDirMatchTimes{i}{j});
            PostUnDirMatchTimes{i}{j} = PostUnDirMatchTimes{i}{j}(SortedIndices);
            PostUnDirMatchLabels{i}{j} = PostUnDirMatchLabels{i}{j}(SortedIndices);
            PostUnDirPlotLabels{i}{j} = PostUnDirPlotLabels{i}{j}(SortedIndices);
        end
    end
end


for i = 1:Parameters.NoPreDays,
    for j = 1:length(Parameters.SyllablePreDirResults{i}{1}),
        PreDirMatchTimes{i}{j} = [];
        PreDirMatchLabels{i}{j} = [];
        PreDirPlotLabels{i}{j} = [];
        for k = 1:length(SyllLabels),
            if (~isempty(Parameters.SyllablePreDirResults{i}{k}{j}))
                Indices = find(Parameters.SyllablePreDirResults{i}{k}{j}(:,1) >= Results.ActualSyllableTemplateMatchThreshold{k});
                PreDirMatchTimes{i}{j} = [PreDirMatchTimes{i}{j}; Parameters.SyllablePreDirResults{i}{k}{j}(Indices,2)];
                PreDirMatchLabels{i}{j} = [PreDirMatchLabels{i}{j}; char(ones(size(Indices))*SyllLabels(k))];
                PreDirPlotLabels{i}{j} = [PreDirPlotLabels{i}{j}; (ones(size(Indices))*k)];
            end
        end
        if (~isempty(PreDirMatchLabels{i}{j}))
            [SortedVals, SortedIndices] = sort(PreDirMatchTimes{i}{j});
            PreDirMatchTimes{i}{j} = PreDirMatchTimes{i}{j}(SortedIndices);
            PreDirMatchLabels{i}{j} = PreDirMatchLabels{i}{j}(SortedIndices);
            PreDirPlotLabels{i}{j} = PreDirPlotLabels{i}{j}(SortedIndices);
        end
    end
end

for i = 1:Parameters.NoPostDays,
    for j = 1:length(Parameters.SyllablePostDirResults{i}{1}),
        PostDirMatchTimes{i}{j} = [];
        PostDirMatchLabels{i}{j} = [];
        PostDirPlotLabels{i}{j} = [];
        for k = 1:length(SyllLabels),
            if (~isempty(Parameters.SyllablePostDirResults{i}{k}{j}))
                Indices = find(Parameters.SyllablePostDirResults{i}{k}{j}(:,1) >= Results.ActualSyllableTemplateMatchThreshold{k});
                PostDirMatchTimes{i}{j} = [PostDirMatchTimes{i}{j}; Parameters.SyllablePostDirResults{i}{k}{j}(Indices,2)];
                PostDirMatchLabels{i}{j} = [PostDirMatchLabels{i}{j}; char(ones(size(Indices))*SyllLabels(k))];
                PostDirPlotLabels{i}{j} = [PostDirPlotLabels{i}{j}; (ones(size(Indices))*k)];
            end
        end
        if (~isempty(PostDirMatchLabels{i}{j}))
            [SortedVals, SortedIndices] = sort(PostDirMatchTimes{i}{j});
            PostDirMatchTimes{i}{j} = PostDirMatchTimes{i}{j}(SortedIndices);
            PostDirMatchLabels{i}{j} = PostDirMatchLabels{i}{j}(SortedIndices);
            PostDirPlotLabels{i}{j} = PostDirPlotLabels{i}{j}(SortedIndices);
        end
    end
end


figure;
Index = 1;
for i = 1:length(PreUnDirMatchTimes{1}),
    Matches = find(PreUnDirMatchLabels{1}{i} == 'd');
    hold on;
    for j = 1:length(Matches),
        scatter(PreUnDirMatchTimes{1}{i} - PreUnDirMatchTimes{1}{i}(Matches(j)), ones(size(PreUnDirMatchTimes{1}{i}))*Index, 25, PreUnDirPlotLabels{1}{i}, 'filled');
        Index = Index + 1;
    end
end
axis tight;
temp = axis;
axis([-1.5 1.5 1 temp(end)]);

figure;
Index = 1;
for i = 1:length(PostUnDirMatchTimes{end}),
    Matches = find(PostUnDirMatchLabels{end}{i} == 'd');
    hold on;
    for j = 1:length(Matches),
        scatter(PostUnDirMatchTimes{end}{i} - PostUnDirMatchTimes{end}{i}(Matches(j)), ones(size(PostUnDirMatchTimes{end}{i}))*Index, 25, PostUnDirPlotLabels{end}{i}, 'filled');
        Index = Index + 1;
    end
end
axis tight;
temp = axis;
axis([-1.5 1.5 1 temp(end)]);

figure;
Index = 1;
for i = 1:length(PreDirMatchTimes{1}),
    Matches = find(PreDirMatchLabels{1}{i} == 'd');
    hold on;
    for j = 1:length(Matches),
        scatter(PreDirMatchTimes{1}{i} - PreDirMatchTimes{1}{i}(Matches(j)), ones(size(PreDirMatchTimes{1}{i}))*Index, 25, PreDirPlotLabels{1}{i}, 'filled');
        Index = Index + 1;
    end
end
axis tight;
temp = axis;
axis([-1.5 1.5 1 temp(end)]);

figure;
Index = 1;
for i = 1:length(PostDirMatchTimes{end}),
    Matches = find(PostDirMatchLabels{end}{i} == 'd');
    hold on;
    for j = 1:length(Matches),
        scatter(PostDirMatchTimes{end}{i} - PostDirMatchTimes{end}{i}(Matches(j)), ones(size(PostDirMatchTimes{end}{i}))*Index, 25, PostDirPlotLabels{end}{i}, 'filled');
        Index = Index + 1;
    end
end
axis tight;
temp = axis;
axis([-1.5 1.5 1 temp(end)]);

disp('Finished Analysis');