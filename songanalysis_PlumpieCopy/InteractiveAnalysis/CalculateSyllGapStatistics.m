function [DirSongLengths, DirSyllables, DirGaps, UnDirSongLengths, UnDirSyllables, UnDirGaps, MedianMotif] = CalculateSyllGapStatistics(DirFileInfo, UnDirFileInfo, MotifString)

Motifs = struct('Onsets', [], 'Offsets', [], 'Labels', [], 'FileName', []);

Index = 1;
DirIndex = 1;

DirSongLengths = [];
DirSyllables = [];
DirGaps = [];

UnDirSongLengths = [];
UnDirSyllables = [];
UnDirGaps = [];

if (~isempty(DirFileInfo))
    for i = 1:length(DirFileInfo.NoteOnsets),
        Matches = strfind(DirFileInfo.NoteLabels{i}, MotifString);
        for j = 1:length(Matches),
            Motifs(Index).Onsets = DirFileInfo.NoteOnsets{i}(Matches(j):(Matches(j) + length(MotifString) - 1));
            Motifs(Index).Offsets = DirFileInfo.NoteOffsets{i}(Matches(j):(Matches(j) + length(MotifString) - 1));
            Motifs(Index).Labels = DirFileInfo.NoteLabels{i}(Matches(j):(Matches(j) + length(MotifString) - 1));
            Motifs(Index).FileName = DirFileInfo.FileNames{i};

            DirSongLengths(DirIndex,:) = Motifs(Index).Offsets(end) - Motifs(Index).Onsets(1);
            DirSyllables.Start(DirIndex,:) = Motifs(Index).Onsets;
            DirSyllables.End(DirIndex,:) = Motifs(Index).Offsets;
            DirSyllables.Length(DirIndex,:) = DirSyllables.End(DirIndex,:) - DirSyllables.Start(DirIndex,:);
            DirSyllables.Index(DirIndex) = i;
            
            DirGaps.Start(DirIndex,:) = Motifs(Index).Offsets(1:(end - 1));
            DirGaps.End(DirIndex,:) = Motifs(Index).Onsets(2:end);
            DirGaps.Length(DirIndex,:) = DirGaps.End(DirIndex,:) - DirGaps.Start(DirIndex,:);
            DirGaps.Index(DirIndex) = i;
            
            DirIndex = DirIndex + 1;
            Index = Index + 1;
        end
    end
    disp(['No of directed song motifs is ', num2str(DirIndex - 1)]);
end

UnDirIndex = 1;
if (~isempty(UnDirFileInfo))
    for i = 1:length(UnDirFileInfo.NoteOnsets),
        Matches = strfind(UnDirFileInfo.NoteLabels{i}, MotifString);
        for j = 1:length(Matches),
            Motifs(Index).Onsets = UnDirFileInfo.NoteOnsets{i}(Matches(j):(Matches(j) + length(MotifString) - 1));
            Motifs(Index).Offsets = UnDirFileInfo.NoteOffsets{i}(Matches(j):(Matches(j) + length(MotifString) - 1));
            Motifs(Index).Labels = UnDirFileInfo.NoteLabels{i}(Matches(j):(Matches(j) + length(MotifString) - 1));
            Motifs(Index).FileName = UnDirFileInfo.FileNames{i};

            UnDirSongLengths(UnDirIndex,:) = Motifs(Index).Offsets(end) - Motifs(Index).Onsets(1);            
            UnDirSyllables.Start(UnDirIndex,:) = Motifs(Index).Onsets;
            UnDirSyllables.End(UnDirIndex,:) = Motifs(Index).Offsets;
            UnDirSyllables.Length(UnDirIndex,:) = UnDirSyllables.End(UnDirIndex,:) - UnDirSyllables.Start(UnDirIndex,:);
            UnDirSyllables.Index(UnDirIndex) = i;
            
            UnDirGaps.Start(UnDirIndex,:) = Motifs(Index).Offsets(1:(end - 1));
            UnDirGaps.End(UnDirIndex,:) = Motifs(Index).Onsets(2:end);
            UnDirGaps.Length(UnDirIndex,:) = UnDirGaps.End(UnDirIndex,:) - UnDirGaps.Start(UnDirIndex,:);
            UnDirGaps.Index(UnDirIndex) = i;            
            
            UnDirIndex = UnDirIndex + 1;
            Index = Index + 1;
        end
    end
    disp(['No of undirected song motifs is ', num2str(UnDirIndex - 1)]);    
end


for i = 1:length(Motifs),
    MotifDurations(i) = Motifs(i).Offsets(end) - Motifs(i).Onsets(1);
end

[SortedDurations, SortedIndices] = sort(MotifDurations);

if (mod(length(SortedDurations), 2) == 0)
    MedianIndex = SortedIndices(length(SortedIndices)/2);
else
    MedianIndex = SortedIndices((length(SortedIndices) + 1)/2);
end

MedianMotif.FileName{1} = Motifs(MedianIndex).FileName;
MedianMotif.StartTime = Motifs(MedianIndex).Onsets(1);
MedianMotif.SyllableLengths = Motifs(MedianIndex).Offsets - Motifs(MedianIndex).Onsets;
MedianMotif.GapLengths = Motifs(MedianIndex).Onsets(2:end) - Motifs(MedianIndex).Offsets(1:(end - 1));
MedianMotif.Length = MotifDurations(MedianIndex);
MedianMotif.SyllableStartings = Motifs(MedianIndex).Onsets - MedianMotif.StartTime;
MedianMotif.GapStartings = Motifs(MedianIndex).Offsets(1:(end - 1)) - MedianMotif.StartTime;

disp('Calculated syllable gap statistics');
