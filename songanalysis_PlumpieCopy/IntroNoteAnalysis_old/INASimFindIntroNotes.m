function [SimMotifIntroNotes] = INASimFindIntroNotes(Labels, MotifSimIndices)

SimMotifIntroNotes = [];

BoutNo = 1;
MotifNo = 1;
Motifs = MotifSimIndices;
for j = 1:length(Motifs),
    if (j == 1)
        IntroNotes = find(Labels(1:Motifs(j)) == 'i');
        if (~isempty(IntroNotes))
            SimMotifIntroNotes.Indices{MotifNo} = IntroNotes;
            SimMotifIntroNotes.NoofINs{MotifNo} = length(IntroNotes);
            MotifNo = MotifNo + 1;
        else
            SimMotifIntroNotes.NoofINs{MotifNo} = 0;
            MotifNo = MotifNo + 1;
        end
    else
        IntroNotes = find(Labels(Motifs(j-1):Motifs(j)) == 'i');
        if (~isempty(IntroNotes))
            IntroNotes = IntroNotes + Motifs(j-1) - 1;
            SimMotifIntroNotes.Indices{MotifNo} = IntroNotes;
            SimMotifIntroNotes.NoofINs{MotifNo} = length(IntroNotes);
            MotifNo = MotifNo + 1;
        else
            SimMotifIntroNotes.NoofINs{MotifNo} = 0;
            MotifNo = MotifNo + 1;
        end
    end
end
