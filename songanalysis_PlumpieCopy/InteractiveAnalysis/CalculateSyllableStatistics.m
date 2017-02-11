function [Syllables, Gaps] = CalculateSyllableStatistics(FileInfo, Motif)

Songs = strfind(FileInfo.NoteLabels,Motif);
if (isempty(Songs))
    return;
end

Syllables = cell(length(Songs),1);
Gaps = cell(length(Songs), 1);
    
for SongNo = 1:length(Songs),
    for Syllable = 1:length(Motif),
        Syllables{SongNo}.Start{Syllable} = FileInfo.NoteOnsets{SongNo}(Songs{SongNo} + (Syllable - 1));
        Syllables{SongNo}.End{Syllable} = FileInfo.NoteOffsets{SongNo}(Songs{SongNo} + (Syllable - 1));
        Syllables{SongNo}.Length{Syllable} = Syllables{SongNo}.End{Syllable} - Syllables{SongNo}.Start{Syllable};

        if (Syllable ~= length(Motif))
            Gaps{SongNo}.Start{Syllable} = FileInfo.NoteOffsets{SongNo}(Songs{SongNo} + (Syllable - 1));
            Gaps{SongNo}.End{Syllable} = FileInfo.NoteOnsets{SongNo}(Songs{SongNo} + (Syllable));
            Gaps{SongNo}.Length{Syllable} = Gaps{SongNo}.End{Syllable} - Gaps{SongNo}.Start{Syllable};
        end
    end
end