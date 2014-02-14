function [Raster, PST, SpikeTrain, SpikeWaveforms] = WarpSpikeTrains(FileInfo, Motif, MedianMotif, BinSize, Latency)

Raster = [];
PST = [];
SpikeTrain = [];
SpikeWaveforms = [];

PreDataTime = 0.05;

Edges = -PreDataTime:BinSize:(MedianMotif.Length);
Songs = strfind(FileInfo.Notes.NoteLabels,Motif);
if (length(Songs) == 0)
    return;
end

for SongNo = 1:length(Songs),
    TempSpikes = [];
    TempSpikeAmplitudes = [];
    TempSpikeWaveforms = [];
    OnsetSpikes = FileInfo.SpikeData.Times(find((FileInfo.SpikeData.Times > (FileInfo.Syllables.Start(SongNo,1) - PreDataTime)) & (FileInfo.SpikeData.Times < FileInfo.Syllables.Start(SongNo,1))));
    OnsetSpikeAmplitudes = FileInfo.SpikeData.Amplitudes(find((FileInfo.SpikeData.Times > (FileInfo.Syllables.Start(SongNo,1) - PreDataTime)) & (FileInfo.SpikeData.Times < FileInfo.Syllables.Start(SongNo,1))));
    OnsetSpikeWaveforms = FileInfo.SpikeData.Waveforms(find((FileInfo.SpikeData.Times > (FileInfo.Syllables.Start(SongNo,1) - PreDataTime)) & (FileInfo.SpikeData.Times < FileInfo.Syllables.Start(SongNo,1))),:);    
    if (length(OnsetSpikes > 0))
        TempSpikes = [TempSpikes; (OnsetSpikes - FileInfo.Syllables.Start(SongNo,1))];
        TempSpikeAmplitudes = [TempSpikeAmplitudes; OnsetSpikeAmplitudes];
        TempSpikeWaveforms = [TempSpikeWaveforms; OnsetSpikeWaveforms];
    end
    
    for Syllable = 1:length(Motif),
        SyllableSpikes = FileInfo.SpikeData.Times(find((FileInfo.SpikeData.Times > FileInfo.Syllables.Start(SongNo,Syllable)) & (FileInfo.SpikeData.Times < FileInfo.Syllables.End(SongNo,Syllable))));
        SyllableSpikeAmplitudes = FileInfo.SpikeData.Amplitudes(find((FileInfo.SpikeData.Times > FileInfo.Syllables.Start(SongNo,Syllable)) & (FileInfo.SpikeData.Times < FileInfo.Syllables.End(SongNo,Syllable))));
        SyllableSpikeWaveforms = FileInfo.SpikeData.Waveforms(find((FileInfo.SpikeData.Times > FileInfo.Syllables.Start(SongNo,Syllable)) & (FileInfo.SpikeData.Times < FileInfo.Syllables.End(SongNo,Syllable))),:);        
        if (length(SyllableSpikes > 0))
            SyllableSpikes = (((SyllableSpikes - FileInfo.Syllables.Start(SongNo,Syllable)) * MedianMotif.SyllableLengths(Syllable))/(FileInfo.Syllables.Length(SongNo, Syllable))) + MedianMotif.SyllableStartings(Syllable);
            TempSpikes = [TempSpikes; SyllableSpikes];
            TempSpikeAmplitudes = [TempSpikeAmplitudes; SyllableSpikeAmplitudes];
            TempSpikeWaveforms = [TempSpikeWaveforms; SyllableSpikeWaveforms];
        end

        if (Syllable ~=length(Motif))
            GapSpikes = FileInfo.SpikeData.Times(find((FileInfo.SpikeData.Times > FileInfo.Gaps.Start(SongNo,Syllable)) & (FileInfo.SpikeData.Times < FileInfo.Gaps.End(SongNo,Syllable))));
            GapSpikeAmplitudes = FileInfo.SpikeData.Amplitudes(find((FileInfo.SpikeData.Times > FileInfo.Gaps.Start(SongNo,Syllable)) & (FileInfo.SpikeData.Times < FileInfo.Gaps.End(SongNo,Syllable))));            
            GapSpikeWaveforms = FileInfo.SpikeData.Waveforms(find((FileInfo.SpikeData.Times > FileInfo.Gaps.Start(SongNo,Syllable)) & (FileInfo.SpikeData.Times < FileInfo.Gaps.End(SongNo,Syllable))),:);                        
            if (length(GapSpikes > 0))
                GapSpikes = (((GapSpikes - FileInfo.Gaps.Start(SongNo,Syllable)) * MedianMotif.GapLengths(Syllable))/(FileInfo.Gaps.Length(SongNo, Syllable))) + MedianMotif.GapStartings(Syllable);
                TempSpikes = [TempSpikes; GapSpikes];
                TempSpikeAmplitudes = [TempSpikeAmplitudes; GapSpikeAmplitudes];                
                TempSpikeWaveforms = [TempSpikeWaveforms; GapSpikeWaveforms];
            end
        end
    end
    
%    OffsetSpikes = FileInfo.SpikeData.Times(find((FileInfo.SpikeData.Times > FileInfo.Syllables.End(SongNo,end)) & (FileInfo.SpikeData.Times < (FileInfo.Syllables.End(SongNo,end) + 0.2))));
%    OffsetSpikeAmplitudes = FileInfo.SpikeData.Amplitudes(find((FileInfo.SpikeData.Times > FileInfo.Syllables.End(SongNo,end)) & (FileInfo.SpikeData.Times < (FileInfo.Syllables.End(SongNo,end) + 0.2))));    
%    OffsetSpikeWaveforms = FileInfo.SpikeData.Waveforms(find((FileInfo.SpikeData.Times > FileInfo.Syllables.End(SongNo,end)) & (FileInfo.SpikeData.Times < (FileInfo.Syllables.End(SongNo,end) + 0.2))),:);        
%    if (length(OffsetSpikes > 0))
%        OffsetSpikes = OffsetSpikes - FileInfo.Syllables.End(SongNo,end) + MedianMotif.Length;
%        TempSpikes = [TempSpikes; OffsetSpikes];
%        TempSpikeAmplitudes = [TempSpikeAmplitudes; OffsetSpikeAmplitudes];                        
%        TempSpikeWaveforms = [TempSpikeWaveforms; OffsetSpikeWaveforms];
%    end
    
    if (length(TempSpikes) > 0)
        Raster = [Raster; [TempSpikes (ones(length(TempSpikes),1) * (SongNo-1)) TempSpikeAmplitudes]];
        SpikeTrain{SongNo} = TempSpikes(find((TempSpikes > (-Latency)) & (TempSpikes < MedianMotif.Length))) + Latency;
        SpikeWaveforms = [SpikeWaveforms; TempSpikeWaveforms];
    end
    if (length(TempSpikes) > 0)
        PST(SongNo,:) = histc(TempSpikes,Edges);
    else
        PST(SongNo,:) = histc([100 100], Edges);
    end
end

PST = PST/BinSize;
