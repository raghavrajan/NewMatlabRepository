function [] = CompareSyllableAmplitudeMethods(DirectoryName, FileList, NoteFileDir, FileType, MotifSylls, Baseline)
 
PresentDir = pwd;
cd(DirectoryName);

% First get all the files
Fid = fopen(FileList, 'r');
Files = textscan(Fid, '%s', 'DeLimiter', '\n');
Files = Files{1};
fclose(Fid);

% Now load up the note files
for i = 1:length(Files),
    NoteInfo{i} = load(fullfile(DirectoryName, NoteFileDir, [Files{i}, '.not.mat']));
end

% Now for each of the files, calculate log amplitude in all the different
% ways and then calculate peak, peak-to-peak and mean values

for i = 1:length(MotifSylls),
    SyllIndex(i) = 0;
end

FHigh = 860;
FLow = 8600;
FFTWinSize = 5; % in ms

fprintf('Total of %d files\n', length(Files));
fprintf('\n');
for i = 1:length(Files),
    fprintf('%d > ', i);
    [RawData, Fs] = GetData(DirectoryName, Files{i}, FileType, 0);
    FiltData = bandpass(RawData, Fs, FHigh, FLow);
    
    LogAmplitudeKao = ASSLCalculateLogAmplitudeKao(RawData, Fs, [], [], [], FHigh, FLow);
    LogAmplitudeAronovFee = ASSLCalculateLogAmplitudeAronovFee(RawData, Fs, [], FFTWinSize, [], FHigh, FLow);
    
    [m_spec_deriv , m_AM, m_FM ,m_Entropy , m_amplitude ,m_Freq, m_PitchGoodness , m_Pitch , Pitch_chose , Pitch_weight] = deriv(FiltData, Fs);
    LogAmplitudeSAP = m_amplitude;
    
    % First Raw Waveform
    for j = 1:length(MotifSylls),
        Matches = find(NoteInfo{i}.labels == MotifSylls(j));
        for MatchNum = Matches(:)',
            SyllIndex(j) = SyllIndex(j) + 1;
            StartIndex = round((NoteInfo{i}.onsets(MatchNum)/1000) * Fs);
            EndIndex = round((NoteInfo{i}.offsets(MatchNum)/1000) * Fs);
            if (StartIndex < 1)
                StartIndex = 1;
            end
            if (EndIndex > length(FiltData))
                EndIndex = length(FiltData);
            end
            Sylls{j}.MeanRMSAmplitude(SyllIndex(j)) = mean(sqrt(sum(FiltData(StartIndex:EndIndex).^2)));
            Sylls{j}.PeakAmplitude(SyllIndex(j)) = max(FiltData(StartIndex:EndIndex));
            Sylls{j}.PeaktoPeakAmplitude(SyllIndex(j)) = max(FiltData(StartIndex:EndIndex)) - min(FiltData(StartIndex:EndIndex));
            
            Sylls{j}.MeanKaoAmplitude(SyllIndex(j)) = mean(LogAmplitudeKao(StartIndex:EndIndex));
            Sylls{j}.PeakKaoAmplitude(SyllIndex(j)) = max(LogAmplitudeKao(StartIndex:EndIndex));
            
            Sylls{j}.MeanAronovFeeAmplitude(SyllIndex(j)) = mean(LogAmplitudeAronovFee(StartIndex:EndIndex));
            Sylls{j}.PeakAronovFeeAmplitude(SyllIndex(j)) = max(LogAmplitudeAronovFee(StartIndex:EndIndex));
            
            % Now for SAP, I have to adjust the time base as the time base
            % is different for SAP amplitude
            SAPTime = linspace(1/Fs, length(RawData)/Fs, length(LogAmplitudeSAP));
            SAPFs = 1/(SAPTime(2) - SAPTime(1));
            StartIndex = round((NoteInfo{i}.onsets(MatchNum)/1000) * SAPFs);
            EndIndex = round((NoteInfo{i}.offsets(MatchNum)/1000) * SAPFs);
            if (StartIndex < 1)
                StartIndex = 1;
            end
            if (EndIndex > length(LogAmplitudeSAP))
                EndIndex = length(LogAmplitudeSAP);
            end
            
            Sylls{j}.MeanSAPAmplitude(SyllIndex(j)) = mean(LogAmplitudeSAP(StartIndex:EndIndex));
            Sylls{j}.PeakSAPAmplitude(SyllIndex(j)) = max(LogAmplitudeSAP(StartIndex:EndIndex));
        end
    end
end

% Now plot the figures for correlation between all of these
FieldNames = {'MeanRMSAmplitude' 'MeanKaoAmplitude' 'MeanSAPAmplitude' 'MeanAronovFeeAmplitude' 'PeakSAPAmplitude' 'PeakAronovFeeAmplitude' 'PeakAmplitude' 'PeaktoPeakAmplitude' 'PeakKaoAmplitude'};
for SyllNum = 1:length(MotifSylls),
    for i = 1:length(FieldNames),
        for j = 1:length(FieldNames),
            [r{SyllNum}(i,j), p{SyllNum}(i,j)] = corr(eval(['[Sylls{', num2str(SyllNum), '}.', FieldNames{i},'(:)]']), eval(['[Sylls{', num2str(SyllNum), '}.', FieldNames{j},'(:)]'])); 
        end
    end
    figure;
    r{SyllNum}(find(p{SyllNum} >= 0.05)) = NaN;
    imagesc(r{SyllNum});
    set(gca, 'FontSize', 14);
    set(gca, 'YTick', 1:1:length(FieldNames), 'YTickLabel', FieldNames);
    set(gca, 'XTick', 1:1:length(FieldNames), 'XTickLabel', FieldNames);
    set(gca, 'XTickLabelRotation', 45);
    title([FileList, ' : Syll ', MotifSylls(SyllNum)]);
    set(gcf, 'Color', 'w');
    set(gcf, 'Position', [680 300 900 700]);
    colorbar
end

% Now just for mean rms amplitude and mean aronov fee log amplitude, I have
% to figure out the effect of adjusting boundaries. I will adjust
% boundaries in a couple of ways:
% 1. Use a common threshold based on checking all waveforms
% 2. Use peaks/troughs in the derivative

fprintf('\n');
for i = 1:length(MotifSylls),
    SyllRawAmp{i} = [];
    SyllLogAmplitude{i} = [];
    SyllRawAmpMeanValues{i} = [];
    SyllLogAmplitudeMeanValues{i} = [];
    ActualStartIndex{i} = [];
    ActualEndIndex{i} = [];
end

PrePostPadding = 0.01; % in seconds

for i = 1:length(Files),
    fprintf('%d > ', i);
    [RawData, Fs] = GetData(DirectoryName, Files{i}, FileType, 0);
    FiltData = bandpass(RawData, Fs, FHigh, FLow);
    
    LogAmplitude = ASSLCalculateLogAmplitudeAronovFee(RawData, Fs, [], FFTWinSize, [], FHigh, FLow);
    LogAmplitude = LogAmplitude - Baseline;
    
    % First Raw Waveform
    for j = 1:length(MotifSylls),
        Matches = find(NoteInfo{i}.labels == MotifSylls(j));
        for MatchNum = Matches(:)',
            SyllIndex(j) = SyllIndex(j) + 1;
            StartIndex = round(((NoteInfo{i}.onsets(MatchNum)/1000) * Fs) - (PrePostPadding * Fs));
            EndIndex = round(((NoteInfo{i}.offsets(MatchNum)/1000) * Fs) + (PrePostPadding * Fs));
            if (StartIndex < 1)
                StartIndex = 1;
                continue;
            end
            if (EndIndex > length(FiltData))
                continue;
                EndIndex = length(FiltData);
            end
            SyllRawAmp{j}{end+1} = FiltData(StartIndex:EndIndex);
            SyllLogAmplitude{j}{end+1} = LogAmplitude(StartIndex:EndIndex);
            
            % Now get the values for these syllables for the raw segmented
            % boundaries
            ActualStartIndex{j}(end+1) = StartIndex;
            ActualEndIndex{j}(end+1) = EndIndex;
            SyllRawAmpMeanValues{j}(end+1,1) = mean(sqrt(sum(FiltData(StartIndex:EndIndex).^2)));
            SyllLogAmplitudeMeanValues{j}(end+1,1) = mean(LogAmplitude(StartIndex:EndIndex));
        end
    end
end

GaussianLen = 3;
Width = 0.005; % in ms
XGauss = 1:1:(1 + round(2 * GaussianLen * Width * (Fs)));
XGauss = XGauss - (length(XGauss) + 1)/2;
GaussWin = (1/((Width * Fs) * sqrt(2 * pi))) * exp(-(XGauss.*XGauss)/(2 * (Width * Fs) * (Width * Fs)));

AllowedAdjust = 0.008; % in seconds
PrePostPaddingIndex = round(PrePostPadding * Fs);
% Now to modify boundaries using either a common threshold or the
% derivative - this can only be done for log amplitude values
for i = 1:length(SyllLogAmplitude),
   PreMinValue = [];
   PostMinValue = [];
   for j = 1:length(SyllLogAmplitude{i}),
       PreMinValue(j) = min(SyllLogAmplitude{i}{j}(1:PrePostPaddingIndex));
       PostMinValue(j) = min(SyllLogAmplitude{i}{j}(end-PrePostPaddingIndex+1:end));
   end
   PreCommonThreshold = max(PreMinValue);
   PostCommonThreshold = max(PostMinValue);
   % Now adjust boundaries according to new threshold and calculate mean
   % values
   for j = 1:length(SyllLogAmplitude{i}),
       if (SyllLogAmplitude{i}{j}(PrePostPaddingIndex + 1) > PreCommonThreshold)
           [Val, NewStartIndex] = find(SyllLogAmplitude{i}{j}(1:PrePostPaddingIndex) <= PreCommonThreshold, 1, 'last');
       else
           [Val, NewStartIndex] = find(SyllLogAmplitude{i}{j}(PrePostPaddingIndex+1:end) >= PreCommonThreshold, 1, 'first');
           if (~isempty(NewStartIndex))
               NewStartIndex = NewStartIndex + PrePostPaddingIndex;
           end
       end
       if isempty(NewStartIndex)
           NewStartIndex = PrePostPaddingIndex + 1;
       end
       
       if (SyllLogAmplitude{i}{j}(end - PrePostPaddingIndex) > PostCommonThreshold)
           [Val, NewEndIndex] = find(SyllLogAmplitude{i}{j}(end - PrePostPaddingIndex + 1:end) <= PostCommonThreshold, 1, 'first');
           if (~isempty(NewEndIndex))
               NewEndIndex = NewEndIndex + length(SyllLogAmplitude{i}{j}) - PrePostPaddingIndex;
           end
       else
           [Val, NewEndIndex] = find(SyllLogAmplitude{i}{j}(1:end - PrePostPaddingIndex) >= PostCommonThreshold, 1, 'last');
       end
       if isempty(NewEndIndex)
           NewEndIndex = length(SyllLogAmplitude{i}{j}) - PrePostPaddingIndex;
       end
      SyllLogAmplitudeMeanValues{i}(j,2) = mean(SyllLogAmplitude{i}{j}(NewStartIndex:NewEndIndex));
      SyllRawAmpMeanValues{i}(j,2) = mean(sqrt(sum(SyllRawAmp{i}{j}(NewStartIndex:NewEndIndex).^2)));
      
      % Now to find peaks and troughs in the derivative and use that also
      % to adjust boundaries. I will also use the condition that the peaks
      % and troughs have to be > (mean + 1.5*std) or < (mean - 1.5*std)
      % respectively

      ActualStartIndex{i}(j) = PrePostPaddingIndex + 1;
      ActualEndIndex{i}(j) = length(SyllLogAmplitude{i}{j}) - PrePostPaddingIndex;

      LogAmpDerivative = diff(conv(SyllLogAmplitude{i}{j}, GaussWin, 'same'));
      %Threshold = [(mean(LogAmpDerivative) + 1.5*std(LogAmpDerivative)) (mean(LogAmpDerivative) - 1.5*std(LogAmpDerivative))];
      Threshold = prctile(LogAmpDerivative, [5 95]);
      [PosPks, PosPkIndices] = findpeaks(LogAmpDerivative,'MinPeakHeight', Threshold(1)); 
      % find closes PosPk next to actual boundary - allow for only 5ms
      % adjustment maximum
      if (isempty(PosPkIndices))
          PosPkIndex = NaN;
          AdjustPosValue{i}(j) = NaN;
      else
          [MinVal, MinPkIndex] = min(abs(PosPkIndices - ActualStartIndex{i}(j)));
          AdjustPosValue{i}(j) = MinVal;
          if (abs(MinVal) <= round(AllowedAdjust * Fs))
              PosPkIndex = PosPkIndices(MinPkIndex);
          else
              PosPkIndex = NaN;
          end
      end
      
      [NegPks, NegPkIndices] = findpeaks(-LogAmpDerivative, 'MinPeakHeight', -Threshold(2));
      % find closes NegPk next to actual boundary - allow for only 5ms
      % adjustment maximum
      if (isempty(NegPkIndices))
          NegPkIndex = NaN;
          AdjustNegValue{i}(j) = NaN;
      else
          [MinVal, MinPkIndex] = min(abs(NegPkIndices - ActualEndIndex{i}(j)));
          AdjustNegValue{i}(j) = MinVal;
          if (abs(MinVal) <= round(AllowedAdjust * Fs))
              NegPkIndex = NegPkIndices(MinPkIndex);
          else
              NegPkIndex = NaN;
          end
      end
      
      if (~isnan(NegPkIndex) && ~isnan(PosPkIndex))
          SyllLogAmplitudeMeanValues{i}(j,3) = mean(SyllLogAmplitude{i}{j}(PosPkIndex:NegPkIndex));
          SyllRawAmpMeanValues{i}(j,3) = mean(sqrt(sum(SyllRawAmp{i}{j}(PosPkIndex:NegPkIndex).^2)));
      else
          SyllLogAmplitudeMeanValues{i}(j,3) = NaN;
          SyllRawAmpMeanValues{i}(j,3) = NaN;
      end
      
      if (j <= 6)
          figure;
          hold on;
          set(gcf, 'Position', [83 558 1834 420]);
          [Ax, H1, H2] = plotyy(1:1:length(SyllLogAmplitude{i}{j}), SyllLogAmplitude{i}{j}, 1:1:length(SyllLogAmplitude{i}{j})-1, LogAmpDerivative);
          axes(Ax(1));
          set(H1, 'MarkerSize', 4, 'Color', 'b', 'Marker', 'o');
          axis tight;
          Temp = axis;
          plot([Temp(1) PrePostPaddingIndex], ones(1,2)*PreCommonThreshold, 'k--');
          text(mean([Temp(1) PrePostPaddingIndex]), 1.05*PreCommonThreshold, num2str(PreCommonThreshold));

          plot([(Temp(2) - PrePostPaddingIndex) Temp(2)], ones(1,2)*PostCommonThreshold, 'k--');
          text(mean([(Temp(2) - PrePostPaddingIndex) Temp(2)]), 1.05*PostCommonThreshold, num2str(PostCommonThreshold));

          plot(ones(1,2)*PrePostPaddingIndex+1, Temp(3:4), 'k--');
          plot(ones(1,2)*(length(SyllLogAmplitude{i}{j}) - PrePostPaddingIndex), Temp(3:4), 'k--');
          
          plot(NewStartIndex, SyllLogAmplitude{i}{j}(NewStartIndex), 'rs', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
          plot(NewEndIndex, SyllLogAmplitude{i}{j}(NewEndIndex), 'rs', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
          
          if (~isnan(NegPkIndex) && ~isnan(PosPkIndex))
            plot(PosPkIndex, SyllLogAmplitude{i}{j}(PosPkIndex), 'cs', 'MarkerSize', 8, 'MarkerFaceColor', 'c');
            plot(NegPkIndex, SyllLogAmplitude{i}{j}(NegPkIndex), 'cs', 'MarkerSize', 8, 'MarkerFaceColor', 'c');
          end
          axes(Ax(2));
          set(H2, 'MarkerSize', 4, 'Color', 'm', 'Marker', 'o');
          axis tight;
      end
   end
end

% Display % of syllables that are NaNs for the derivative method
fprintf('\n');
disp(['% of syllables for which boundaries could not be adjusted with a common threshold']);
for i = 1:length(SyllLogAmplitudeMeanValues),
    disp([num2str(i), ': ', num2str(length(find(isnan(SyllLogAmplitudeMeanValues{i}(:,2)))) * 100 / size(SyllLogAmplitudeMeanValues{i},1)), ' %']);
end
disp(['% of syllables for which boundaries could not be adjusted with a derivative']);
for i = 1:length(SyllLogAmplitudeMeanValues),
    disp([num2str(i), ': ', num2str(length(find(isnan(SyllLogAmplitudeMeanValues{i}(:,3)))) * 100 / size(SyllLogAmplitudeMeanValues{i},1)), ' %']);
end

% Display the adjusted values mean and range
for i = 1:length(SyllLogAmplitudeMeanValues),
    disp([num2str(i), ': Positive adjust: median - ', num2str(nanmedian(AdjustPosValue{i})/Fs), ': iqr - ', num2str(iqr(AdjustPosValue{i})/Fs), ' s']);
    disp([num2str(i), ': Negative adjust: median - ', num2str(nanmedian(AdjustNegValue{i})/Fs), ': iqr - ', num2str(iqr(AdjustNegValue{i})/Fs), ' s']);
end

% Correlation between different measures of amplitude
for i = 1:length(SyllLogAmplitudeMeanValues),
    disp(['Syllable #', num2str(i)]);
    disp(['Correlation between unadjusted and adjusted with common threshold = ', num2str(corr(SyllLogAmplitudeMeanValues{i}(:,1), SyllLogAmplitudeMeanValues{i}(:,2), 'rows', 'complete'))]);
    disp(['Correlation between unadjusted and adjusted with first derivative = ', num2str(corr(SyllLogAmplitudeMeanValues{i}(:,3), SyllLogAmplitudeMeanValues{i}(:,2), 'rows', 'complete'))]);
    disp(['Correlation between adjusted with first derivative and adjusted with common threshold = ', num2str(corr(SyllLogAmplitudeMeanValues{i}(:,1), SyllLogAmplitudeMeanValues{i}(:,3), 'rows', 'complete'))]);
end
disp('Finished comparing syllable amplitudes');