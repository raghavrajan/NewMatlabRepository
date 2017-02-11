function [MedianMotif] = CalculateMedianMotifLengthOneCondition(FileInfo,Motif)

MedianMotif = [];
SongLengths = [];

Songs = strfind(FileInfo.Notes.NoteLabels,Motif);
for i = 1:length(Songs),
    SongLengths(i,:) = FileInfo.Notes.NoteOffsets(Songs(i) + (length(Motif) - 1)) - FileInfo.Notes.NoteOnsets(Songs(i));
end

[SortedSongLengths, SortedSongIndices] = sort(SongLengths);

if (mod(length(SongLengths),2) == 0)
    MedianMotifIndex = SortedSongIndices(length(SortedSongIndices)/2);
else
    MedianMotifIndex = SortedSongIndices((length(SortedSongIndices) + 1)/2);
end

FileIndex = find(cumsum(FileInfo.RecordLengths) < FileInfo.Syllables.Start(MedianMotifIndex,1),1,'last');
if (length(FileIndex) == 0)
    FileIndex = 0;
end

MedianMotif.FileName = FileInfo.FileNames(FileIndex + 1);
MedianMotif.StartTime = FileInfo.Syllables.Start(MedianMotifIndex,1) - sum(FileInfo.RecordLengths(1:FileIndex));

for i = 1:length(Motif),
    MedianMotif.SyllableLengths(i) = FileInfo.Notes.NoteOffsets(Songs(MedianMotifIndex) + (i-1)) - FileInfo.Notes.NoteOnsets(Songs(MedianMotifIndex) + (i-1));
    if (i ~= length(Motif))
        MedianMotif.GapLengths(i) = FileInfo.Notes.NoteOnsets(Songs(MedianMotifIndex) + (i)) - FileInfo.Notes.NoteOffsets(Songs(MedianMotifIndex) + (i-1));
    end
end

MedianMotif.Length = SongLengths(MedianMotifIndex);

MedianMotif.SyllableStartings(1) = 0;

for i = 2:length(Motif),
    MedianMotif.SyllableStartings(i) = MedianMotif.SyllableStartings(i-1) + MedianMotif.SyllableLengths(i-1) + MedianMotif.GapLengths(i-1);
end

for i = 1:(length(Motif)-1),
    MedianMotif.GapStartings(i) = MedianMotif.SyllableStartings(i) + MedianMotif.SyllableLengths(i);
end
