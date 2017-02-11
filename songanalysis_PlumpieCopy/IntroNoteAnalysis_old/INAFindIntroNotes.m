function [MotifIntroNotes, NonMotifIntroNotes] = INAFindIntroNotes(BoutIndices, MotifIndices, Onsets, Offsets, Labels, InterBoutInt, FileLength)

MotifIntroNotes = [];
NonMotifIntroNotes = [];
NonMotifIntroNotes.Indices = [];

BoutNo = 0;
MotifNo = 1;
for i = 1:2:length(BoutIndices),
    TempIntroNotes = [];
    ActualIntroNotes = [];
    BoutNo = BoutNo + 1;
    MotifIndex = find((MotifIndices >= BoutIndices(i)) & (MotifIndices <= BoutIndices(i+1)));
    if (isempty(MotifIndex))
        continue;
    else
        Motifs = MotifIndices(MotifIndex);
        if (Onsets(Motifs(1)) <= InterBoutInt)
            continue;
        else
            for j = 1:length(Motifs),
                MotifIntroNotes.MotifOnsets{MotifNo} = Onsets(Motifs(j));
                if (j == 1)
                    IntroNotes = find(Labels(BoutIndices(i):Motifs(j)) == 'i');
                    if (~isempty(IntroNotes))
                        IntroNotes = IntroNotes + BoutIndices(i) - 1;
                        Intervals = Onsets(IntroNotes(2:end)) - Offsets(IntroNotes(1:end-1));
                        LongIntervals = find(Intervals > InterBoutInt);
                        if (~isempty(LongIntervals))
                            IntroNotes = IntroNotes((LongIntervals(end) + 1):end);
                        end
                        IntroNotes = sort(IntroNotes);
                        INIndexIntervals = diff(IntroNotes);
                        LongIntervals = find(INIndexIntervals > 1);
                        if (~isempty(LongIntervals))
                            IntroNotes = IntroNotes((LongIntervals(end) + 1):end);
                        end
                            
                        ActualIntroNotes = [ActualIntroNotes IntroNotes];
                        MotifIntroNotes.Indices{MotifNo} = IntroNotes;
                        MotifIntroNotes.NoofINs{MotifNo} = length(IntroNotes);
                        MotifIntroNotes.Durs{MotifNo} = Offsets(IntroNotes) - Onsets(IntroNotes);
                        MotifIntroNotes.Intervals{MotifNo} = Onsets(IntroNotes + 1) - Offsets(IntroNotes);
                        MotifIntroNotes.Onsets{MotifNo} = Onsets(IntroNotes);
                        MotifIntroNotes.Offsets{MotifNo} = Offsets(IntroNotes);
                        if (IntroNotes(1) == 1)
                            MotifIntroNotes.VocalInterval{MotifNo} = Onsets(IntroNotes(1));
                        else
                            MotifIntroNotes.VocalInterval{MotifNo} = Onsets(IntroNotes(1)) - Offsets(IntroNotes(1) - 1);
                        end
                        MotifIntroNotes.BoutNo{MotifNo} = BoutNo;
                        MotifNo = MotifNo + 1;
                    else
                        MotifIntroNotes.NoofINs{MotifNo} = 0;
                        if (Motifs(j) == 1)
                            MotifIntroNotes.VocalInterval{MotifNo} = Onsets(Motifs(j));
                        else
                            MotifIntroNotes.VocalInterval{MotifNo} = Onsets(Motifs(j)) - Offsets(Motifs(j) - 1);
                        end
                        MotifIntroNotes.BoutNo{MotifNo} = BoutNo;
                        MotifNo = MotifNo + 1;
                    end
                else
                    IntroNotes = find(Labels(Motifs(j-1):Motifs(j)) == 'i');
                    if (~isempty(IntroNotes))
                        IntroNotes = IntroNotes + Motifs(j-1) - 1;
                        Intervals = Onsets(IntroNotes(2:end)) - Offsets(IntroNotes(1:end-1));
                        LongIntervals = find(Intervals > InterBoutInt);
                        if (~isempty(LongIntervals))
                            IntroNotes = IntroNotes((LongIntervals(end) + 1):end);
                        end
                        IntroNotes = sort(IntroNotes);
                        INIndexIntervals = diff(IntroNotes);
                        LongIntervals = find(INIndexIntervals > 1);
                        if (~isempty(LongIntervals))
                            IntroNotes = IntroNotes((LongIntervals(end) + 1):end);
                        end
                        ActualIntroNotes = [ActualIntroNotes IntroNotes];
                        MotifIntroNotes.Indices{MotifNo} = IntroNotes;
                        MotifIntroNotes.NoofINs{MotifNo} = length(IntroNotes);
                        MotifIntroNotes.Durs{MotifNo} = Offsets(IntroNotes) - Onsets(IntroNotes);
                        MotifIntroNotes.Intervals{MotifNo} = Onsets(IntroNotes + 1) - Offsets(IntroNotes);
                        MotifIntroNotes.Onsets{MotifNo} = Onsets(IntroNotes);
                        MotifIntroNotes.Offsets{MotifNo} = Offsets(IntroNotes);
                        MotifIntroNotes.VocalInterval{MotifNo} = Onsets(IntroNotes(1)) - Offsets(IntroNotes(1) - 1);
                        MotifIntroNotes.BoutNo{MotifNo} = BoutNo;
                        MotifNo = MotifNo + 1;
                    else
                        MotifIntroNotes.NoofINs{MotifNo} = 0;
                        MotifIntroNotes.VocalInterval{MotifNo} = Onsets(Motifs(j)) - Offsets(Motifs(j) - 1);
                        MotifIntroNotes.BoutNo{MotifNo} = BoutNo;
                        MotifNo = MotifNo + 1;
                    end
                end
            end
            AllIntroNotes = find(Labels(BoutIndices(i):BoutIndices(i+1)) == 'i');
            AllIntroNotes = AllIntroNotes + BoutIndices(i) - 1;
            NonMotifIntroNotes.Indices = [NonMotifIntroNotes.Indices setdiff(AllIntroNotes, ActualIntroNotes)];
        end
    end
end
NonMotifIntroNotes.Durs{1} = Offsets(NonMotifIntroNotes.Indices) - Onsets(NonMotifIntroNotes.Indices);
for i = 1:length(NonMotifIntroNotes.Indices),
    if (NonMotifIntroNotes.Indices(i) == length(Onsets))
        NonMotifIntroNotes.Intervals{1}(i,1) = FileLength - Offsets(NonMotifIntroNotes.Indices(i));
    else
        NonMotifIntroNotes.Intervals{1}(i,1) = Onsets(NonMotifIntroNotes.Indices(i) + 1) - Offsets(NonMotifIntroNotes.Indices(i));
    end
end
NonMotifIntroNotes.Onsets{1} = Onsets(NonMotifIntroNotes.Indices);
NonMotifIntroNotes.Offsets{1} = Offsets(NonMotifIntroNotes.Indices);