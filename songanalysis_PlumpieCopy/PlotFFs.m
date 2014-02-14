function [] = PlotFFs(DirectoryName, BirdName, Syll)

cd(DirectoryName);

Files = dir([BirdName, '*', Syll, '*mat']);

for i = 1:length(Files),
    load(Files(i).name);
    plot((ones(size(cell2mat(ffreq(:,2))))*i), cell2mat(ffreq(:,2)), 'k+');
    hold on;
end

axis tight;
temp = axis;
axis([0 4 temp(3) temp(4)]);
