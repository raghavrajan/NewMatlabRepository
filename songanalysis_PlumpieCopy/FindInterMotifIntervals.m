function [InterMotifIntervals] = FindInterMotifIntervals(DirectoryName, BirdName, Motif);

cd(DirectoryName);

NotesFiles = dir([BirdName, '*.not.mat']);

InterMotifIntervals = [];
for i = 1:length(NotesFiles),
    clear Notes;
    Notes = load(NotesFiles(i).name);
    disp(Notes.labels);
    temp = strfind(Notes.labels, Motif);
    for j = 2:length(temp),
        InterMotifIntervals(end + 1) = Notes.onsets(temp(j)) - Notes.offsets(temp(j) - 1);
    end
end

figure
bar([0:5:500], histc(InterMotifIntervals, [0:5:500]), 'histc');
disp(['Minimum inter motif interval is ', num2str(min(InterMotifIntervals)), ' and the median is ', num2str(median(InterMotifIntervals))]);
