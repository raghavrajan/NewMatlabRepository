function [] = PlotInterBoutInterval_vs_BoutNumber(DirectoryName, FileList, FileType, NoteFileDir, MotifSyllables)

% ========================== Help text ====================================
% This function can be used to screen through labelled files to figure out
% how many bouts there are for a range of inter-bout intervals. The plot
% generated can then be used to choose an optimum inter-bout interval for
% which there is sufficient number of bouts. This optimum inter-bout
% interval can then be entered and the program will generate a text file
% with the list of files that have bouts that satisfy the inter-bout
% interval criterion. These can then be used for further analysis.
% Inputs: 
%   1) DirectoryName - directory where all data files are stored
%   2) FileList - text file with the list of files that have to be checked
%   3) FileType - type of data file
%   4) NoteFileDir - directory where note files can be found
%   5) MotifSyllables - syllables that should be considered as motif
%   syllables. The program will only give all bouts with motif syllables in
%   them
% Assumption is all these files are triggered files. For continuous files
% there should not be any problem as all the data is there and needs to be
% split up later.
% =========================================================================

% Set values for Inter-bout interval
InterboutInterval = 500:250:3000; % in ms

% First open the filelist and get all filenames
Fid = fopen(fullfile(DirectoryName, FileList), 'r');
Files = textscan(Fid, '%s', 'DeLimiter', '\n');
Files = Files{1};
fclose(Fid);

% Now load up the notefile corresponding to each file. Tag on an onset at 0
% and an offset at the file length of the file. Then for each of the
% different inter-bout intervals, find all bouts in the file that have
% motif syllables and also have inter-bout interval >= specified inter-bout
% interval. Right at the beginning if the file does not have any of the
% motif syllables, then the file can be skipped.
for i = 1:length(Files),
    Notes = load(fullfile(NoteFileDir, [Files{i}, '.not.mat']));
    [RawData, Fs] = GetData(DirectoryName, Files{i}, FileType, 0);
    FileLen = length(RawData)/Fs;
    
    % Skip file if it does not have any syllables
    if (isempty(Notes.labels))
        Bouts(i,1:length(InterboutInterval)) = 0;
        continue;
    end
    
    % Skip file if it does not have any motif syllables
    HasMotifSyllables = 0;
    for j = 1:length(MotifSyllables),
        if (~isempty(find(Notes.labels == MotifSyllables(j))))
            HasMotifSyllables = 1;
            break;
        end
    end
    
    if (HasMotifSyllables == 0)
        Bouts(i,1:length(InterboutInterval)) = 0;
        continue;
    end
    
    % Now assume there is an offset at 0
    Notes.offsets = [0; Notes.offsets(:)];
    % and an onset at the end of the file
    Notes.onsets = [Notes.onsets(:); FileLen*1000];
    
    Intervals = Notes.onsets - Notes.offsets;
    
    for j = 1:length(InterboutInterval),
        LongIntervals = find(Intervals >= InterboutInterval(j));
        NumIntervals = length(LongIntervals);
        
        if (NumIntervals < 2)
            Bouts(i,j) = 0;
        else
            BoutNum = 0;
            HasMotifSyllables = 0;
            for k = 1:NumIntervals-1,
                BoutLabels = Notes.labels(LongIntervals(k):(LongIntervals(k+1)-1));
                for MotifSyllIndex = 1:length(MotifSyllables),
                        if (~isempty(find(BoutLabels == MotifSyllables(MotifSyllIndex))))
                            HasMotifSyllables = 1;
                            break;
                        end
                end
                if (HasMotifSyllables == 1)
                    BoutNum = BoutNum + 1;
                end
            end
            Bouts(i,j) = BoutNum;
        end
    end
end

% Now plot the number of bouts vs. Inter-bout interval
figure;
set(gcf, 'Color', 'w');
plot(InterboutInterval, sum(Bouts), 'ko-');
for i = 1:length(InterboutInterval),
    text(InterboutInterval(i), sum(Bouts(:,i))*1.05, ['(', num2str(sum(Bouts(:,i))), ')'], 'FontSize', 8);
end
xlabel('Inter-bout interval (ms)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Num bouts', 'FontSize', 14, 'FontWeight', 'bold');

% Now choose required inter-bout interval based on number of bouts needed
InterboutIntervalChoice = inputdlg('Choose the inter-bout interval that you want in ms', 'Inter-bout interval choice');
InterboutIntervalChoice = str2double(InterboutIntervalChoice{1});

Index = find(InterboutInterval == InterboutIntervalChoice);
ValidFiles = find(Bouts(:,Index));
Fid = fopen(fullfile(DirectoryName, [FileList, '.', num2str(InterboutInterval(Index)), 'ms.txt']), 'w');
for i = 1:length(ValidFiles),
    fprintf(Fid, '%s\n', Files{ValidFiles(i)});
end
fclose(Fid);

disp('Finished');