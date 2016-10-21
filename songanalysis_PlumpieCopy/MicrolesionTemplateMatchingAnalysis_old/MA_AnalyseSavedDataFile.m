function [Results, Params] = MA_AnalyseSavedDataFile(Parameters, PaddingLength)

%====== First analyse bout length for pre and post ========================
% Bout length is stored for each day in an array
% Row indicates a day and the first three columns are # of bouts, mean bout
% length and standard deviation of bout length for dir song. The next 3
% columns are the same for undir song. 
% First set of rows are for pre days and the next set of rows are for post
% days

for i = 1:Parameters.NoPreDays,
    DirBouts = cell2mat(Parameters.PreDirBoutLens{i});
    UnDirBouts = cell2mat(Parameters.PreUnDirBoutLens{i});
    
    Results.BoutLength.Dir.NumBouts(i) = length(DirBouts);
    Results.BoutLength.Dir.Mean(i) = mean(DirBouts);
    Results.BoutLength.Dir.STD(i) = std(DirBouts);
    Results.BoutLength.Dir.Median(i) = median(DirBouts);
    Results.BoutLength.Dir.IQR(i) = iqr(DirBouts);
    
    Results.BoutLength.UnDir.NumBouts(i) = length(UnDirBouts);
    Results.BoutLength.UnDir.Mean(i) = mean(UnDirBouts);
    Results.BoutLength.UnDir.STD(i) = std(UnDirBouts);
    Results.BoutLength.UnDir.Median(i) = median(UnDirBouts);
    Results.BoutLength.UnDir.IQR(i) = iqr(UnDirBouts);

    clear DirBouts UnDirBouts;
end

for i = 1:Parameters.NoPostDays,
    DirBouts = cell2mat(Parameters.PostDirBoutLens{i});
    UnDirBouts = cell2mat(Parameters.PostUnDirBoutLens{i});
    
    Results.BoutLength.Dir.NumBouts(i + Parameters.NoPreDays) = length(DirBouts);
    Results.BoutLength.Dir.Mean(i + Parameters.NoPreDays) = mean(DirBouts);
    Results.BoutLength.Dir.STD(i + Parameters.NoPreDays) = std(DirBouts);
    Results.BoutLength.Dir.Median(i + Parameters.NoPreDays) = median(DirBouts);
    Results.BoutLength.Dir.IQR(i + Parameters.NoPreDays) = iqr(DirBouts);
    
    Results.BoutLength.UnDir.NumBouts(i + Parameters.NoPreDays) = length(UnDirBouts);
    Results.BoutLength.UnDir.Mean(i + Parameters.NoPreDays) = mean(UnDirBouts);
    Results.BoutLength.UnDir.STD(i + Parameters.NoPreDays) = std(UnDirBouts);
    Results.BoutLength.UnDir.Median(i + Parameters.NoPreDays) = median(UnDirBouts);
    Results.BoutLength.UnDir.IQR(i + Parameters.NoPreDays) = iqr(UnDirBouts);
    
    clear DirBouts UnDirBouts;
end

%==========================================================================

%====== Analyse rhythm for each bout of singing ===========================
% For each bout, take the syllable onsets and offsets and convert each bout
% into a series of 1s and 0s corresponding to the onset and offset period
% respectively - sampling rate of 1kHz. This can then be used for 
% generating the rhythm spectrogram

% It is easy to do syllable and gap duration histograms as part of this for
% loop, so am going to include that too in this.


Params.DiscreteBout_Fs = 500; % sampling rate for representation of bouts as 1s and 0s - Hz
Params.NoofBins = 50;
Params.DurEdges = linspace(-2.5, 0, Params.NoofBins);
Params.NFFT = 2^nextpow2(round(2*Params.DiscreteBout_Fs));
Params.Freq = Params.DiscreteBout_Fs/2*linspace(0, 1, Params.NFFT/2+1);
Params.FreqIndices = find(Params.Freq <= 30);
% Params.HanningWinLen = 0.01; % in sec
% Params.HanningWin = hanning(round(Params.HanningWinLen * Parameters.PreDirFs{1}{1}));
Params.HanningWin = hanning(Params.NFFT);
Params.HanningWin = Params.HanningWin/sum(Params.HanningWin);

for i = 1:Parameters.NoPreDays,
    DirSyllDurs = [];
    DirGapDurs = [];
    
    DirSyllSyllDurs = [];
    DirSyllGapDurs = [];
    DirGapSyllDurs = [];
    DirGapGapDurs = [];
    
    BoutIndex = 1;
    for j = 1:length(Parameters.PreDirBoutLens{i}),
        [RawData, Fs] = ASSLGetRawData(Parameters.PreDataDir{i}, Parameters.PreDirSongFileNames{i}{j}, Parameters.FileType, 0);
        for k = 1:length(Parameters.PreDirBoutLens{i}{j}),
            SyllableIndices = find((Parameters.PreDirOnsets{i}{j} >= Parameters.PreDirBoutOnsets{i}{j}(k)) & (Parameters.PreDirOnsets{i}{j} <= Parameters.PreDirBoutOffsets{i}{j}(k)));
            
            Results.NumSyllablesPerBout.Dir{i}(BoutIndex) = length(SyllableIndices);
            
            SyllDurs = Parameters.PreDirOffsets{i}{j}(SyllableIndices) - Parameters.PreDirOnsets{i}{j}(SyllableIndices);
            GapDurs = Parameters.PreDirOnsets{i}{j}(SyllableIndices(2:end)) - Parameters.PreDirOffsets{i}{j}(SyllableIndices(1:end-1));
            SyllDurs = SyllDurs(:);
            GapDurs = GapDurs(:);
            
            DirSyllDurs = [DirSyllDurs; SyllDurs(:)];
            DirGapDurs = [DirGapDurs; GapDurs(:)];
            
            DirSyllSyllDurs = [DirSyllSyllDurs; [SyllDurs(1:end-1) SyllDurs(2:end)]]; 
            DirSyllGapDurs = [DirSyllGapDurs; [SyllDurs(1:end-1) GapDurs]]; 
            DirGapSyllDurs = [DirGapSyllDurs; [GapDurs SyllDurs(2:end)]]; 
            DirGapGapDurs = [DirGapGapDurs; [GapDurs(1:end-1) GapDurs(2:end)]]; 
            
            clear SyllDurs GapDurs;
            
            SongBout = RawData(ceil(Parameters.PreDirBoutOnsets{i}{j}(k)*Fs/1000):floor(Parameters.PreDirBoutOffsets{i}{j}(k)*Fs/1000));
            
            % Do SAP Feature calculations
            BoutSyllOnsets = (Parameters.PreDirOnsets{i}{j}(SyllableIndices) - Parameters.PreDirBoutOnsets{i}{j}(k))/1000;
            BoutSyllOffsets = (Parameters.PreDirOffsets{i}{j}(SyllableIndices) - Parameters.PreDirBoutOnsets{i}{j}(k))/1000;
            [Results.Feats.Dir{i}{BoutIndex}, RawFeats] = ASSLCalculateSAPFeatsWithOnsets(SongBout, (1:1:length(SongBout))/Fs, Fs, BoutSyllOnsets, BoutSyllOffsets);
            
            % Do Rhythm spectrogram calculations
            SongBout = SongBout(floor(PaddingLength*Fs):(length(SongBout) - ceil(PaddingLength * Fs)));
            [TempBoutFFT] = MA_CalcRhythmSpectrogram(SongBout, Fs, Params.DiscreteBout_Fs, Params.Freq(Params.FreqIndices)); 
            DirBoutFFT{i}(BoutIndex,:) = TempBoutFFT;
            
            clear TempBoutFFT;
            
            DirOnsetsOffsets{BoutIndex} = zeros(round(Parameters.PreDirBoutLens{i}{j}(k)/1000 * Params.DiscreteBout_Fs), 1);
            for SyllNo = 1:length(SyllableIndices),
                SyllOnset = ceil((Parameters.PreDirOnsets{i}{j}(SyllableIndices(SyllNo)) - Parameters.PreDirBoutOnsets{i}{j}(k))/1000 * Params.DiscreteBout_Fs);
                SyllOffset = floor((Parameters.PreDirOffsets{i}{j}(SyllableIndices(SyllNo)) - Parameters.PreDirBoutOnsets{i}{j}(k))/1000 * Params.DiscreteBout_Fs);
                DirOnsetsOffsets{BoutIndex}(SyllOnset:SyllOffset) = 1;
            end
            BoutIndex = BoutIndex + 1;
        end
    end
    
    UnDirSyllDurs = [];
    UnDirGapDurs = [];
    
    UnDirSyllSyllDurs = [];
    UnDirSyllGapDurs = [];
    UnDirGapSyllDurs = [];
    UnDirGapGapDurs = [];
    
    BoutIndex = 1;
    for j = 1:length(Parameters.PreUnDirBoutLens{i}),
        [RawData, Fs] = ASSLGetRawData(Parameters.PreDataDir{i}, Parameters.PreUnDirSongFileNames{i}{j}, Parameters.FileType, 0);
        for k = 1:length(Parameters.PreUnDirBoutLens{i}{j}),
            SyllableIndices = find((Parameters.PreUnDirOnsets{i}{j} >= Parameters.PreUnDirBoutOnsets{i}{j}(k)) & (Parameters.PreUnDirOnsets{i}{j} <= Parameters.PreUnDirBoutOffsets{i}{j}(k)));
            
            Results.NumSyllablesPerBout.UnDir{i}(BoutIndex) = length(SyllableIndices);
            
            SyllDurs = Parameters.PreUnDirOffsets{i}{j}(SyllableIndices) - Parameters.PreUnDirOnsets{i}{j}(SyllableIndices);
            GapDurs = Parameters.PreUnDirOnsets{i}{j}(SyllableIndices(2:end)) - Parameters.PreUnDirOffsets{i}{j}(SyllableIndices(1:end-1));
            UnDirSyllDurs = [UnDirSyllDurs; SyllDurs(:)];
            UnDirGapDurs = [UnDirGapDurs; GapDurs(:)];
            
            UnDirSyllSyllDurs = [UnDirSyllSyllDurs; [SyllDurs(1:end-1) SyllDurs(2:end)]]; 
            UnDirSyllGapDurs = [UnDirSyllGapDurs; [SyllDurs(1:end-1) GapDurs]]; 
            UnDirGapSyllDurs = [UnDirGapSyllDurs; [GapDurs SyllDurs(2:end)]]; 
            UnDirGapGapDurs = [UnDirGapGapDurs; [GapDurs(1:end-1) GapDurs(2:end)]]; 
            
            clear SyllDurs GapDurs;
            
            SongBout = RawData(ceil(Parameters.PreUnDirBoutOnsets{i}{j}(k)*Fs/1000):floor(Parameters.PreUnDirBoutOffsets{i}{j}(k)*Fs/1000));
            SongBout = SongBout(floor(PaddingLength*Fs):(length(SongBout) - ceil(PaddingLength * Fs)));
            [TempBoutFFT] = MA_CalcRhythmSpectrogram(SongBout, Fs, Params.DiscreteBout_Fs, Params.Freq(Params.FreqIndices)); 
            UnDirBoutFFT{i}(BoutIndex,:) = TempBoutFFT;
            
            clear TempBoutFFT;
            
            UnDirOnsetsOffsets{BoutIndex} = zeros(round(Parameters.PreUnDirBoutLens{i}{j}(k)/1000 * Params.DiscreteBout_Fs), 1);
            for SyllNo = 1:length(SyllableIndices),
                SyllOnset = ceil((Parameters.PreUnDirOnsets{i}{j}(SyllableIndices(SyllNo)) - Parameters.PreUnDirBoutOnsets{i}{j}(k))/1000 * Params.DiscreteBout_Fs);
                SyllOffset = floor((Parameters.PreUnDirOffsets{i}{j}(SyllableIndices(SyllNo)) - Parameters.PreUnDirBoutOnsets{i}{j}(k))/1000 * Params.DiscreteBout_Fs);
                UnDirOnsetsOffsets{BoutIndex}(SyllOnset:SyllOffset) = 1;
            end
            BoutIndex = BoutIndex + 1;
        end
    end
    
    Results.DirSyllDurs{i} = DirSyllDurs;
    Results.DirGapDurs{i} = DirGapDurs;
    Results.UnDirSyllDurs{i} = UnDirSyllDurs;
    Results.UnDirGapDurs{i} = UnDirGapDurs;
    
    DirSyllDurs = log10(DirSyllDurs/1000);
    DirGapDurs = log10(DirGapDurs/1000);
    UnDirSyllDurs = log10(UnDirSyllDurs/1000);
    UnDirGapDurs = log10(UnDirGapDurs/1000);
    
    Results.DirSyllDurHist(i,:) = histc(DirSyllDurs, Params.DurEdges)/length(DirSyllDurs);
    Results.UnDirSyllDurHist(i,:) = histc(UnDirSyllDurs, Params.DurEdges)/length(UnDirSyllDurs);
    Results.DirGapDurHist(i,:) = histc(DirGapDurs, Params.DurEdges)/length(DirGapDurs);
    Results.UnDirGapDurHist(i,:) = histc(UnDirGapDurs, Params.DurEdges)/length(UnDirGapDurs);    
end


for i = 1:Parameters.NoPostDays,
    DirSyllDurs = [];
    DirGapDurs = [];
    BoutIndex = 1;
    for j = 1:length(Parameters.PostDirBoutLens{i}),
        [RawData, Fs] = ASSLGetRawData(Parameters.PostDataDir{i}, Parameters.PostDirSongFileNames{i}{j}, Parameters.FileType, 0);
        for k = 1:length(Parameters.PostDirBoutLens{i}{j}),
            SyllableIndices = find((Parameters.PostDirOnsets{i}{j} >= Parameters.PostDirBoutOnsets{i}{j}(k)) & (Parameters.PostDirOnsets{i}{j} <= Parameters.PostDirBoutOffsets{i}{j}(k)));
            
            Results.NumSyllablesPerBout.Dir{i + Parameters.NoPreDays}(BoutIndex) = length(SyllableIndices);
            
            SyllDurs = Parameters.PostDirOffsets{i}{j}(SyllableIndices) - Parameters.PostDirOnsets{i}{j}(SyllableIndices);
            GapDurs = Parameters.PostDirOnsets{i}{j}(SyllableIndices(2:end)) - Parameters.PostDirOffsets{i}{j}(SyllableIndices(1:end-1));
            DirSyllDurs = [DirSyllDurs; SyllDurs(:)];
            DirGapDurs = [DirGapDurs; GapDurs(:)];
            clear SyllDurs GapDurs;
            
            SongBout = RawData(ceil(Parameters.PostDirBoutOnsets{i}{j}(k)*Fs/1000):floor(Parameters.PostDirBoutOffsets{i}{j}(k)*Fs/1000));
            SongBout = SongBout(floor(PaddingLength*Fs):(length(SongBout) - ceil(PaddingLength * Fs)));
            [TempBoutFFT] = MA_CalcRhythmSpectrogram(SongBout, Fs, Params.DiscreteBout_Fs, Params.Freq(Params.FreqIndices));  
            DirBoutFFT{i + Parameters.NoPreDays}(BoutIndex,:) = TempBoutFFT;
            clear Freq TempBoutFFT;
           
            DirOnsetsOffsets{BoutIndex} = zeros(round(Parameters.PostDirBoutLens{i}{j}(k)/1000 * Params.DiscreteBout_Fs), 1);
            for SyllNo = 1:length(SyllableIndices),
                SyllOnset = ceil((Parameters.PostDirOnsets{i}{j}(SyllableIndices(SyllNo)) - Parameters.PostDirBoutOnsets{i}{j}(k))/1000 * Params.DiscreteBout_Fs);
                SyllOffset = floor((Parameters.PostDirOffsets{i}{j}(SyllableIndices(SyllNo)) - Parameters.PostDirBoutOnsets{i}{j}(k))/1000 * Params.DiscreteBout_Fs);
                DirOnsetsOffsets{BoutIndex}(SyllOnset:SyllOffset) = 1;
            end
            BoutIndex = BoutIndex + 1;
        end
    end
    
    UnDirSyllDurs = [];
    UnDirGapDurs = [];
    BoutIndex = 1;
    for j = 1:length(Parameters.PostUnDirBoutLens{i}),
        [RawData, Fs] = ASSLGetRawData(Parameters.PostDataDir{i}, Parameters.PostUnDirSongFileNames{i}{j}, Parameters.FileType, 0);
        for k = 1:length(Parameters.PostUnDirBoutLens{i}{j}),
            SyllableIndices = find((Parameters.PostUnDirOnsets{i}{j} >= Parameters.PostUnDirBoutOnsets{i}{j}(k)) & (Parameters.PostUnDirOnsets{i}{j} <= Parameters.PostUnDirBoutOffsets{i}{j}(k)));
            
            Results.NumSyllablesPerBout.UnDir{i + Parameters.NoPreDays}(BoutIndex) = length(SyllableIndices);
            
            SyllDurs = Parameters.PostUnDirOffsets{i}{j}(SyllableIndices) - Parameters.PostUnDirOnsets{i}{j}(SyllableIndices);
            GapDurs = Parameters.PostUnDirOnsets{i}{j}(SyllableIndices(2:end)) - Parameters.PostUnDirOffsets{i}{j}(SyllableIndices(1:end-1));
            UnDirSyllDurs = [UnDirSyllDurs; SyllDurs(:)];
            UnDirGapDurs = [UnDirGapDurs; GapDurs(:)];
            clear SyllDurs GapDurs;
            
            
            SongBout = RawData(ceil(Parameters.PostUnDirBoutOnsets{i}{j}(k)*Fs/1000):floor(Parameters.PostUnDirBoutOffsets{i}{j}(k)*Fs/1000));
            SongBout = SongBout(floor(PaddingLength*Fs):(length(SongBout) - ceil(PaddingLength * Fs)));
            [TempBoutFFT] = MA_CalcRhythmSpectrogram(SongBout, Fs, Params.DiscreteBout_Fs, Params.Freq(Params.FreqIndices));   
            UnDirBoutFFT{i + Parameters.NoPreDays}(BoutIndex,:) = TempBoutFFT;
            clear TempBoutFFT;
            
            UnDirOnsetsOffsets{BoutIndex} = zeros(round(Parameters.PostUnDirBoutLens{i}{j}(k)/1000 * Params.DiscreteBout_Fs), 1);
            for SyllNo = 1:length(SyllableIndices),
                SyllOnset = ceil((Parameters.PostUnDirOnsets{i}{j}(SyllableIndices(SyllNo)) - Parameters.PostUnDirBoutOnsets{i}{j}(k))/1000 * Params.DiscreteBout_Fs);
                SyllOffset = floor((Parameters.PostUnDirOffsets{i}{j}(SyllableIndices(SyllNo)) - Parameters.PostUnDirBoutOnsets{i}{j}(k))/1000 * Params.DiscreteBout_Fs);
                UnDirOnsetsOffsets{BoutIndex}(SyllOnset:SyllOffset) = 1;
            end
            BoutIndex = BoutIndex + 1;
        end
    end
    
    Results.DirSyllDurs{i + Parameters.NoPreDays} = DirSyllDurs;
    Results.DirGapDurs{i + Parameters.NoPreDays} = DirGapDurs;
    Results.UnDirSyllDurs{i + Parameters.NoPreDays} = UnDirSyllDurs;
    Results.UnDirGapDurs{i + Parameters.NoPreDays} = UnDirGapDurs;
    
    DirSyllDurs = log10(DirSyllDurs/1000);
    DirGapDurs = log10(DirGapDurs/1000);
    UnDirSyllDurs = log10(UnDirSyllDurs/1000);
    UnDirGapDurs = log10(UnDirGapDurs/1000);
    
    Results.DirSyllDurHist(i + Parameters.NoPreDays,:) = histc(DirSyllDurs, Params.DurEdges)/length(DirSyllDurs);
    Results.UnDirSyllDurHist(i + Parameters.NoPreDays,:) = histc(UnDirSyllDurs, Params.DurEdges)/length(UnDirSyllDurs);
    Results.DirGapDurHist(i + Parameters.NoPreDays,:) = histc(DirGapDurs, Params.DurEdges)/length(DirGapDurs);
    Results.UnDirGapDurHist(i + Parameters.NoPreDays,:) = histc(UnDirGapDurs, Params.DurEdges)/length(UnDirGapDurs);    
end

Results.DirBoutFFT = DirBoutFFT;
Results.UnDirBoutFFT = UnDirBoutFFT;

for i = 1:length(Results.NumSyllablesPerBout.Dir),
    Results.NumSyllablesPerBout.Mean.Dir(i) = mean(Results.NumSyllablesPerBout.Dir{i});
    Results.NumSyllablesPerBout.STD.Dir(i) = std(Results.NumSyllablesPerBout.Dir{i});
    Results.NumSyllablesPerBout.Median.Dir(i) = median(Results.NumSyllablesPerBout.Dir{i});
    Results.NumSyllablesPerBout.IQR.Dir(i) = iqr(Results.NumSyllablesPerBout.Dir{i});
    
    Results.NumSyllablesPerBout.Mean.UnDir(i) = mean(Results.NumSyllablesPerBout.UnDir{i});
    Results.NumSyllablesPerBout.STD.UnDir(i) = std(Results.NumSyllablesPerBout.UnDir{i});
    Results.NumSyllablesPerBout.Median.UnDir(i) = median(Results.NumSyllablesPerBout.UnDir{i});
    Results.NumSyllablesPerBout.IQR.UnDir(i) = iqr(Results.NumSyllablesPerBout.UnDir{i});
end

% Now to calculate different parameters related to the rhythm spectrogram
% based on Goldberg and Fee (2011) and Kaie and Hahnloser (2011)

Index = find(Params.Freq(Params.FreqIndices) >= 3, 1, 'first');

for i = 1:length(Results.DirBoutFFT),
    TempBoutFFT = mean(DirBoutFFT{i});
    LinearFit_Coeffs = polyfit(Params.Freq(Params.FreqIndices(Index:end)), TempBoutFFT(Index:end), 1);
    LinearFit = polyval(LinearFit_Coeffs, Params.Freq(Params.FreqIndices(Index:end)));
    Results.DirRhythm_Modulation(i) = sqrt(sum((LinearFit - TempBoutFFT(Index:end)).^2));
    Results.DirRhythmPeak(i) = max(TempBoutFFT(Index:end));
    
    TempBoutFFT = mean(UnDirBoutFFT{i});
    LinearFit_Coeffs = polyfit(Params.Freq(Params.FreqIndices(Index:end)), TempBoutFFT(Index:end), 1);
    LinearFit = polyval(LinearFit_Coeffs, Params.Freq(Params.FreqIndices(Index:end)));
    Results.UnDirRhythm_Modulation(i) = sqrt(sum((LinearFit - TempBoutFFT(Index:end)).^2));
    Results.UnDirRhythmPeak(i) = max(TempBoutFFT(Index:end));
end

% Now to store the entropy-based measure of variability (Goldberg and Fee,
% JNeurophys, 2011) in a matrix - each row corresponds to a day and there
% are two values in each row, the first one corresponding to dir and the
% second one to undir. Same for gaps.



for i = 1:size(Results.DirSyllDurHist,1),
    Temp = Results.DirSyllDurHist(i,:).*log10(Results.DirSyllDurHist(i,:));
    Results.SyllVar(i,1) = -sum(Temp(~isnan(Temp)))/log10(Params.NoofBins);
    
    Temp = Results.UnDirSyllDurHist(i,:).*log10(Results.UnDirSyllDurHist(i,:));
    Results.SyllVar(i,2) = -sum(Temp(~isnan(Temp)))/log10(Params.NoofBins);
    
    Temp = Results.DirGapDurHist(i,:).*log10(Results.DirGapDurHist(i,:));
    Results.GapVar(i,1) = -sum(Temp(~isnan(Temp)))/log10(Params.NoofBins);
    
    Temp = Results.UnDirGapDurHist(i,:).*log10(Results.UnDirGapDurHist(i,:));
    Results.GapVar(i,2) = -sum(Temp(~isnan(Temp)))/log10(Params.NoofBins);
end
disp('Finished');