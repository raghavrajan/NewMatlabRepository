function [] = PlotDifferentTemplateMatchVersions(BirdName)

figure;
hold on;
cd /home/raghav/MicrolesionAnalysisResults_MT_Templates/ShuffledSongComparisons/

Folders = dir([BirdName, '*']);
cd(Folders(1).name);
Files = dir([BirdName, '*']);
MotifIndex = strfind(Files(1).name, '.Motif');
RootFileName = Files(1).name(1:MotifIndex-1);
load(Files(1).name);

plot(Bout.T(1:length(Bout.MaxBoutSeqMatch)), Bout.MaxBoutSeqMatch, 'b');

cd /home/raghav/MicrolesionAnalysisResults_MT_Templates/PartShuffledSongComparisons/

Folders = dir([BirdName, '*']);

Colors = 'rcy';
for i = 1:min([3 length(Folders)]),
    cd(Folders(i).name);
    Files = dir([RootFileName, '*']);
    load(Files(1).name);

    plot(Bout.T(1:length(Bout.MaxBoutSeqMatch)), Bout.MaxBoutSeqMatch, Colors(i));
    cd ..;
end

cd /home/raghav/MicrolesionAnalysisResults_MT_Templates

Folders = dir([BirdName, '*']);

for i = 1:length(Folders),
    cd(Folders(i).name);
    Files = dir([RootFileName, '*']);
    load(Files(1).name);

    plot(Bout.T(1:length(Bout.MaxBoutSeqMatch)), Bout.MaxBoutSeqMatch, 'g');
    cd ..;
end

axis tight;

set(gca, 'Color', 'k'); set(gcf, 'Color', 'k'); set(gca, 'XColor', 'w'); set(gca, 'YColor', 'w'); set(gca, 'FontSize', 16);
xlabel('Time (sec)', 'FontSize', 16);
ylabel('Match value', 'FontSize', 16);
Legend = legend('Completely shuffled', '25% shuffled', '50% shuffled', '75% shuffled', 'normal');
set(Legend, 'Color', 'w');

