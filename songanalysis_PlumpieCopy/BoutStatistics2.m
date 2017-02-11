function [NoofIntroNotes, NoMotifs, INoteDurs, INoteGapDurs, SpikingActivity] = BoutStatistics2(NoteDir, RawDataDir, MainMotif, ExtraMotif, SpikeFileDir, SpikeChan, SpikeClusterNo)

cd(NoteDir);
NoteFiles = dir(['*.not.mat']);
% for i = 1:length(NoteFiles),
%     NoteFiles(i).name = [NoteFiles(i).name(7:end-3), 'not.mat'];
% end

BinSize = 0.005;

PreDur = 150;
PostDur = 150;
Edges = -PreDur:BinSize*1000:PostDur;

InterBoutInterval = 1000;
MinBoutLen = 1000;
INoteDurs = zeros(10000, 15);
INoteGapDurs = zeros(10000, 15);

for i = 1:20,
    SpikingActivity(i).Raster = [];
    SpikingActivity(i).SyllOFFTime = [];
    SpikingActivity(i).NextSyllONTime = [];
    SpikingActivity(i).PrevSyllOFFTime = [];
    SpikingActivity(i).PST = [];
end

for i = 1:length(NoteFiles),
    [rawdata, Fs] = ReadOKrankData(RawDataDir, NoteFiles(i).name(1:end-8), 1);
    FileLength(i) = length(rawdata)/Fs;
    Notes{i} = load(NoteFiles(i).name);
    SpikeTimes{i} = load([SpikeFileDir, '/Chan', num2str(SpikeChan), '_', NoteFiles(i).name(1:end-8), '.spk']);
    SpikeTimes{i} = SpikeTimes{i}(find(SpikeTimes{i}(:,1) == SpikeClusterNo), 2);
    SpikeTimes{i} = SpikeTimes{i}*1000;
end

BoutNo = 1;

for i = 1:length(Notes),
    disp(NoteFiles(i).name);
    disp(Notes{i}.labels);
    Bouts = find(diff(Notes{i}.onsets) > InterBoutInterval);
    if (length(Bouts) > 0)
        Bout = [];
        if (Notes{i}.onsets(1) > InterBoutInterval)
            if ((Notes{i}.offsets(Bouts(1)) - Notes{i}.onsets(1)) > MinBoutLen)
                Bout(1,:) = [1 Bouts(1)];
            end
        end
        if (length(Bouts) > 1)
            for j = 1:length(Bouts)-1,
                if ((Notes{i}.offsets(Bouts(j+1)) - Notes{i}.onsets(Bouts(j) + 1)) > MinBoutLen)
                    Bout(end+1,:) = [(Bouts(j)+1) Bouts(j+1)];
                end
            end
        end
        if (Notes{i}.offsets(end) < (FileLength(i)*1000 - InterBoutInterval))
            if ((Notes{i}.offsets(end) - Notes{i}.onsets(Bouts(end)+1)) > MinBoutLen)
                Bout(end+1,:) = [(Bouts(end) + 1) length(Notes{i}.onsets)];
            end
        end
    else
        Bout = [];
        if ((Notes{i}.onsets(1) > InterBoutInterval) && (Notes{i}.offsets(end) < (FileLength(i)*1000 - InterBoutInterval)) && ((Notes{i}.offsets(end) - Notes{i}.onsets(1)) > MinBoutLen))
            Bout = [1 length(Notes{i}.onsets)];
        end
    end     
    for j = 1:size(Bout,1),
        if (length(strfind(Notes{i}.labels(Bout(j,1):Bout(j,2)), MainMotif)) < 1)
            continue;
        end
        Motifs = [Bout(j,1) (strfind(Notes{i}.labels(Bout(j,1):Bout(j,2)), MainMotif) + Bout(j,1) - 1) Bout(j,2)];
        for k = 2:length(Motifs)-1,
            IntroNotes = find(Notes{i}.labels(Motifs(k-1):Motifs(k)) == 'i');
            NoofIntroNotes(BoutNo) = length(IntroNotes);
            if (length(IntroNotes) > 0)
                IntroNoteDurations = Notes{i}.offsets(Motifs(k-1)+IntroNotes - 1) - Notes{i}.onsets(Motifs(k-1)+IntroNotes-1);
                IntroNoteDurations = IntroNoteDurations(length(IntroNoteDurations):-1:1);
                GapDurations = Notes{i}.onsets(Motifs(k-1)+IntroNotes) - Notes{i}.offsets(Motifs(k-1)+IntroNotes-1);
                GapDurations = GapDurations(length(GapDurations):-1:1);
                INoteDurs(BoutNo,1:length(IntroNoteDurations)) = IntroNoteDurations;
                INoteGapDurs(BoutNo,1:length(GapDurations)) = GapDurations;
                clear IntroNoteDurations GapDurations;
                INoteIndex = 1;
                for INoteNo = length(IntroNotes):-1:1,
                    IntroNoteIndex = Motifs(k-1) + IntroNotes(INoteNo) - 1;
                    INoteTime = Notes{i}.onsets(IntroNoteIndex);
                    INoteSpikeTimes = SpikeTimes{i}(find((SpikeTimes{i} >= (INoteTime - PreDur)) & (SpikeTimes{i} <= (INoteTime + PostDur))));
                    INoteSpikeTimes = INoteSpikeTimes - INoteTime;
                    SpikingActivity(INoteIndex).PST = [SpikingActivity(INoteIndex).PST; histc(INoteSpikeTimes, Edges)'/BinSize];
                    SpikingActivity(INoteIndex).Raster = [SpikingActivity(INoteIndex).Raster; [INoteSpikeTimes ones(size(INoteSpikeTimes))*BoutNo]];
                    SpikingActivity(INoteIndex).SyllOFFTime = [SpikingActivity(INoteIndex).SyllOFFTime; [(Notes{i}.offsets(IntroNoteIndex) - INoteTime) BoutNo]];
                    SpikingActivity(INoteIndex).NextSyllONTime = [SpikingActivity(INoteIndex).NextSyllONTime; [(Notes{i}.onsets(IntroNoteIndex+1) - INoteTime) BoutNo]];
                    if (INoteNo ~= 1)
                        SpikingActivity(INoteIndex).PrevSyllOFFTime = [SpikingActivity(INoteIndex).PrevSyllOFFTime; [(Notes{i}.offsets(IntroNoteIndex-1) - INoteTime) BoutNo]];
                    end
                    INoteIndex = INoteIndex + 1;
                end
            end
            NoMotifs(BoutNo) = 1 + length(strfind(Notes{i}.labels(Motifs(k):Motifs(k+1)), ExtraMotif));
            disp(['No of intro notes is ', num2str(NoofIntroNotes(BoutNo)), ' and the no of motifs is ', num2str(NoMotifs(BoutNo))]);
            BoutNo = BoutNo + 1;
        end
    end
end