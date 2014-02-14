function [Syllables, Gaps] = ProcessNoteFiles(FileInfo, Motif)

Syllables = [];
Gaps = [];

Songs = strfind(FileInfo.Notes.NoteLabels,Motif);
if (length(Songs) == 0)
    return;
end

for SongNo = 1:length(Songs),
    for Syllable = 1:length(Motif),
        Syllables.Start(SongNo, Syllable) = FileInfo.Notes.NoteOnsets(Songs(SongNo) + (Syllable - 1));
        Syllables.End(SongNo,Syllable) = FileInfo.Notes.NoteOffsets(Songs(SongNo) + (Syllable - 1));    
        Syllables.Length(SongNo,Syllable) = Syllables.End(SongNo,Syllable) - Syllables.Start(SongNo,Syllable);
        
        if (Syllable ~= length(Motif))
            Gaps.Start(SongNo,Syllable) = FileInfo.Notes.NoteOffsets(Songs(SongNo) + (Syllable - 1));
            Gaps.End(SongNo,Syllable) = FileInfo.Notes.NoteOnsets(Songs(SongNo) + (Syllable));
            Gaps.Length(SongNo,Syllable) = Gaps.End(SongNo,Syllable) - Gaps.Start(SongNo,Syllable);
        end
    end
end
