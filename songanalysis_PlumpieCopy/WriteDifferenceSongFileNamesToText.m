function [] = WriteDifferenceSongFileNamesToText(TextFile1, TextFile2)

%% Help
% =========================================================================
% Read in two text files; one is bigger and has more files listed in it,
% while the other is a subset of the first 1. Write out all the files that
% are in the bigger list but not in the smaller list to a separate file.
% =========================================================================

%% Read files
Fid = fopen(TextFile1, 'r');
TextFiles1 = textscan(Fid, '%s', 'DeLimiter', '\n');
TextFiles1 = TextFiles1{1};
fclose(Fid);

Fid = fopen(TextFile2, 'r');
TextFiles2 = textscan(Fid, '%s', 'DeLimiter', '\n');
TextFiles2 = TextFiles2{1};
fclose(Fid);

%% Now determine which list is bigger and then find all the ones in the bigger one that are not in the smaller one
if (length(TextFiles1) > length(TextFiles2))
    Difference = setdiff(TextFiles1, TextFiles2);
    OutputFileName = [TextFile1, '.NonSongFiles.txt'];
else
    Difference = setdiff(TextFiles2, TextFiles1);
    OutputFileName = [TextFile2, '.NonSongFiles.txt'];
end

Fid = fopen(OutputFileName, 'w');
for i = 1:length(Difference),
    fprintf(Fid, '%s\n', Difference{i});
end
fclose(Fid);
disp('Finished');