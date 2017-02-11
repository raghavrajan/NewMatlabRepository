function [] = CheckTimeDifferenceBetweenFiles(DataDir, FileList, FileType)

% PreTriggerDuration and PostTriggerDuration in ms

Fid = fopen(FileList, 'r');
Files = textscan(Fid, '%s', 'DeLimiter', '\n');
Files = Files{1};
fclose(Fid);

Index = 1;
NegativeTimeDifferences = [];
for i = 1:length(Files)-1,
    % Remove the extension if there is one for the file
    [PathStr, FileName, Ext] = fileparts(Files{i});
    
    if (exist(fullfile(DataDir, 'ASSLNoteFiles', [Files{i}, '.not.mat']), 'file'))
        load(fullfile(DataDir, 'ASSLNoteFiles', [Files{i}, '.not.mat']));
        ValidLabels = find(labels ~= '0');
        if (~isempty(ValidLabels))
            CurrentFileTime = str2double(FileName(end-5:end-4))*3600*1000 + str2double(FileName(end-3:end-2))*60*1000 + str2double(FileName(end-1:end))*1000; % filetime in ms

            [PathStr, NextFileName, Ext] = fileparts(Files{i+1});
            NextFileTime = str2double(NextFileName(end-5:end-4))*3600*1000 + str2double(NextFileName(end-3:end-2))*60*1000 + str2double(NextFileName(end-1:end))*1000; % filetime in ms

            [RawData, Fs] = GetData(DataDir, Files{i}, FileType, 0);
            [LogAmplitude] = ASSLCalculateLogAmplitudeAronovFee(RawData, Fs, [], 5, []); % 5 is FFTWinSize in ms
            
            OnsetIndex = round(onsets(ValidLabels(1))*Fs/1000);
            if (OnsetIndex < 1)
                OnsetIndex = 1;
            end
            OffsetIndex = round(offsets(ValidLabels(1))*Fs/1000);
            if (OffsetIndex > length(LogAmplitude))
                OffsetIndex = length(LogAmplitude);
            end
            
            FirstSyllAmplitude = mean(LogAmplitude(OnsetIndex:OffsetIndex));
            
            CurrentFileDur = 1000*length(RawData)/Fs;

            TimeDiff = NextFileTime - CurrentFileTime - CurrentFileDur;
            
            if (TimeDiff < 0)
                NegativeTimeDifferences(end+1) = TimeDiff;
                disp(['Negative time difference of ', num2str(TimeDiff/1000), ' sec between file ', FileName, ' of duration ', num2str(CurrentFileDur), ' sec and file ', NextFileName]);
                TimeDiff = abs(TimeDiff);
            end

            FileTimeDifferences(Index,:) = [onsets(ValidLabels(1)) TimeDiff FirstSyllAmplitude];
            
            Index = Index + 1;
        end
    end
end
disp(['Proportion of negative time differences = ', num2str(100 * length(NegativeTimeDifferences)/size(FileTimeDifferences,1)), '%']);
disp(['Proportion of negative time differences <= -1 sec = ', num2str(100 * length(find(NegativeTimeDifferences <= -1000))/length(NegativeTimeDifferences)), '%']);

figure; % Figure with first syllable onset vs. time difference between files
plot(FileTimeDifferences(:,1)/1000, FileTimeDifferences(:,2)/1000, 'ko');
hold on;
set(gcf, 'PaperPositionMode', 'auto');
set(gcf, 'Color', 'w');
set(gcf, 'Position', [198 558 1650 400]);
set(gca, 'YScale', 'log');
xlabel('Time of first syllable in the file (sec)', 'FontSize', 14);
ylabel('Time between files (sec)', 'FontSize', 14);

axis tight;
Temp = axis;
Temp(1) = 0;
axis(Temp);

PreTriggerDuration = inputdlg('Enter the pre-trigger duration in sec', 'Pre Trigger duration');
PreTriggerDuration = str2double(PreTriggerDuration{1});
plot(PreTriggerDuration * ones(1,2), Temp(3:4), 'r--', 'LineWidth', 2);
plot(Temp(1:2), PreTriggerDuration * ones(1,2), 'r--', 'LineWidth', 2);

NumSyllsWithOnsetDuringPreTrigger = find(FileTimeDifferences(:,1) <= PreTriggerDuration);

title([FileList, ': proportion of pre-trigger syllables with file difference <= pre trigger duration = ', num2str(100 * length(find(FileTimeDifferences(NumSyllsWithOnsetDuringPreTrigger,2) <= PreTriggerDuration))/length(NumSyllsWithOnsetDuringPreTrigger)), '%'], 'FontSize', 14);
saveas(gcf, fullfile('/home/raghav/LabNoteBookPlots', [FileList, '.FirstSyllTime.vs.FileTimeDifference.png']), 'png');

figure; % Figure with first syllable onset vs. time difference between files with size and color of marker dependent on amplitude of syllable
set(gcf, 'PaperPositionMode', 'auto');
scatter(FileTimeDifferences(:,1)/1000, FileTimeDifferences(:,2)/1000, ceil((FileTimeDifferences(:,3) + 1.001*abs(min(FileTimeDifferences(:,3))))*2), ceil((FileTimeDifferences(:,3) + 1.001*abs(min(FileTimeDifferences(:,3))))*2))
hold on;
set(gcf, 'Color', 'w');
set(gcf, 'Position', [198 558 1650 400]);
set(gca, 'YScale', 'log');
xlabel('Time of first syllable in the file (sec)', 'FontSize', 12);
ylabel('Time between files (sec)', 'FontSize', 12);

axis tight;
Temp = axis;
Temp(1) = 0;
axis(Temp);
colorbar;

plot(PreTriggerDuration * ones(1,2), Temp(3:4), 'r--', 'LineWidth', 2);
plot(Temp(1:2), PreTriggerDuration * ones(1,2), 'r--', 'LineWidth', 2);
title('First syllable time vs. time since last file: symbol size and color reflects first syllable amplitude', 'FontSize', 14);
text(mean(Temp(1:2)), 1.5*PreTriggerDuration, 'Pre trigger duration (sec)', 'FontSize', 14);

% Now find all the first syllables that fall within a window 50ms before
% the pre-trigger duration to 100ms after the pre-trigger duration and
% calculate the amplitude of these syllables - this will approximately be
% the trigger amplitude.

TriggerWindow = [(PreTriggerDuration - 0.05) (PreTriggerDuration + 0.1)];
TriggerSyllables = find((FileTimeDifferences(:,1)/1000 >= TriggerWindow(1)) & (FileTimeDifferences(:,1)/1000 <= TriggerWindow(2)));
TriggerAmplitude = median(FileTimeDifferences(TriggerSyllables,3));

% Now scale and write this trigger amplitude on previous graph in the top
% right corner.
ScaledTriggerAmplitude = ceil((TriggerAmplitude + + 1.001*abs(min(FileTimeDifferences(:,3))))*2);
text(Temp(1) + 0.25*diff(Temp(1:2)), Temp(3) + 0.75*diff(Temp(3:4))/5, ['Trigger amplitude is ', num2str(ScaledTriggerAmplitude), ' dB'], 'FontSize', 10);

% Now find the pre-trigger syllables that are greater than the trigger
% amplitude and see what proportion of them have file time differences less
% than pre-trigger duration
PreTriggerSyllables_GreaterThanTriggerAmplitude = find((FileTimeDifferences(:,1)/1000 <= TriggerWindow(1)) & (FileTimeDifferences(:,3) >= TriggerAmplitude));
text(Temp(1) + 0.25*diff(Temp(1:2)), Temp(3) + 0.75*diff(Temp(3:4))/10, ['Proportion of pre-trigger syllables with amplitudes >= TriggerAmplitude and file time differences < Pre-trigger duration = ', num2str(100 * length(find(FileTimeDifferences(PreTriggerSyllables_GreaterThanTriggerAmplitude,2)/1000 <= PreTriggerDuration))/length(PreTriggerSyllables_GreaterThanTriggerAmplitude)), '%'], 'FontSize', 10);

saveas(gcf, fullfile('/home/raghav/LabNoteBookPlots', [FileList, '.FirstSyllTime.vs.FileTimeDifference.SyllAmplitudeBasedColorSize.png']), 'png');
disp('Finished');

