function [MedianMotif] = CalculateMedianMotifLength(DirFileInfo,UnDirFileInfo,Motif,RasterPlotFigure)

MedianMotif = [];
DirSongLengths = [];
UnDirSongLengths = [];

DirSongs = strfind(DirFileInfo.Notes.NoteLabels,Motif);
for i = 1:length(DirSongs),
    DirSongLengths(i,:) = DirFileInfo.Notes.NoteOffsets(DirSongs(i) + (length(Motif) - 1)) - DirFileInfo.Notes.NoteOnsets(DirSongs(i));
end

UnDirSongs = strfind(UnDirFileInfo.Notes.NoteLabels,Motif);
for i = 1:length(UnDirSongs),
    UnDirSongLengths(i,:) = UnDirFileInfo.Notes.NoteOffsets(UnDirSongs(i) + (length(Motif) - 1)) - UnDirFileInfo.Notes.NoteOnsets(UnDirSongs(i));
end

SongLengths = [DirSongLengths; UnDirSongLengths];
[SortedSongLengths, SortedSongIndices] = sort(SongLengths);

if (mod(length(SongLengths),2) == 0)
    MedianMotifIndex = SortedSongIndices(length(SortedSongIndices)/2);
else
    MedianMotifIndex = SortedSongIndices((length(SortedSongIndices) + 1)/2);
end

if (MedianMotifIndex > length(DirSongLengths))
    
    MedianMotifIndex = MedianMotifIndex - length(DirSongLengths);  
    
    FileIndex = find(cumsum(UnDirFileInfo.RecordLengths) < UnDirFileInfo.Syllables.Start(MedianMotifIndex,1),1,'last');
    if (length(FileIndex) == 0)
        MedianMotif.FileName = UnDirFileInfo.FileNames(1);
        MedianMotif.StartTime = UnDirFileInfo.Syllables.Start(MedianMotifIndex,1);
    else
        MedianMotif.FileName = UnDirFileInfo.FileNames(FileIndex + 1);
        MedianMotif.StartTime = UnDirFileInfo.Syllables.Start(MedianMotifIndex,1) - sum(UnDirFileInfo.RecordLengths(1:FileIndex));
    end
    
    for i = 1:length(Motif),
        MedianMotif.SyllableLengths(i) = UnDirFileInfo.Notes.NoteOffsets(UnDirSongs(MedianMotifIndex) + (i-1)) - UnDirFileInfo.Notes.NoteOnsets(UnDirSongs(MedianMotifIndex) + (i-1));
        if (i ~= length(Motif))
            MedianMotif.GapLengths(i) = UnDirFileInfo.Notes.NoteOnsets(UnDirSongs(MedianMotifIndex) + (i)) - UnDirFileInfo.Notes.NoteOffsets(UnDirSongs(MedianMotifIndex) + (i-1));
        end
    end

    MedianMotif.Length = UnDirSongLengths(MedianMotifIndex);
else
    
    FileIndex = find(cumsum(DirFileInfo.RecordLengths) < DirFileInfo.Syllables.Start(MedianMotifIndex,1),1,'last');
    if (length(FileIndex) == 0)
        FileIndex = 0;
    end
    
    MedianMotif.FileName = DirFileInfo.FileNames(FileIndex + 1);
    MedianMotif.StartTime = DirFileInfo.Syllables.Start(MedianMotifIndex,1) - sum(DirFileInfo.RecordLengths(1:FileIndex));
    
    for i = 1:length(Motif),
        MedianMotif.SyllableLengths(i) = DirFileInfo.Notes.NoteOffsets(DirSongs(MedianMotifIndex) + (i-1)) - DirFileInfo.Notes.NoteOnsets(DirSongs(MedianMotifIndex) + (i-1));
        if (i ~= length(Motif))
            MedianMotif.GapLengths(i) = DirFileInfo.Notes.NoteOnsets(DirSongs(MedianMotifIndex) + (i)) - DirFileInfo.Notes.NoteOffsets(DirSongs(MedianMotifIndex) + (i-1));
        end
    end

    MedianMotif.Length = DirSongLengths(MedianMotifIndex);
end

MedianMotif.SyllableStartings(1) = 0;

for i = 2:length(Motif),
    MedianMotif.SyllableStartings(i) = MedianMotif.SyllableStartings(i-1) + MedianMotif.SyllableLengths(i-1) + MedianMotif.GapLengths(i-1);
end

for i = 1:(length(Motif)-1),
    MedianMotif.GapStartings(i) = MedianMotif.SyllableStartings(i) + MedianMotif.SyllableLengths(i);
end
