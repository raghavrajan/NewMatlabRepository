function [] = MakeStimulusSetForShikhaWithoutNoise(StimulusFileList, NumRepetitions, MinInterval, MaxInterval, OutputVolume, OutputFileName)

%==========================================================================
% Usage: MakeStimuliForShikha(NumStimuli, StimulusFileList, NumRepetitions, MinInterval, MaxInterval, OutputVolume, OutputFileName)
%
% Inputs:
% StimulusFileList - a text file with list of files (all wav) which have 
%                    stimuli.
%                    From each file, only one stimulus will be taken. Each
%                    row in the file will have a filename, and the onset 
%                    and offset times (in sec) of the stimulus within that
%                    file. The last row should have a file with a 0.5 - 1s
%                    duration corresponding to silence.
% NumRepetitions - number of repetitions of each stimulus
% MinInterval - minimum interval between two consecutive stimuli (in sec)
% MaxInterval - maximum interval between two consecutive stimuli (in sec)
% OutputVolume - a number between 0 and 1 that sets the volume of the ouput
%               stimulus
% OutputFileName - string which will be used for naming all output files
%
% Outputs:
% Does not produce any outputs, but instead writes out wave files for
% stimuli and also writes out a text file that lists the stimuli and the
% times at which they're played.
% 
% =========================================================================

% First read stimulus files and get stimuli and times
Fid = fopen(StimulusFileList, 'r');
Temp = textscan(Fid, '%s', 'DeLimiter', '\n');
Temp = Temp{1};
fclose(Fid);

NumStimuli = length(Temp) - 1;

for i = 1:length(Temp),
    Info = textscan(Temp{i}, '%s%f%f', 'DeLimiter', ' ');
    if (i <= NumStimuli)
        StimFileName{i} = Info{1}{1};
        StimOnsets{i} = Info{2};
        StimOffsets{i} = Info{3};
        [TempData, Fs] = wavread(StimFileName{i});
        Stimulus{i} = TempData(round(Fs*StimOnsets{i}):round(Fs*StimOffsets{i}));
        Stimulus{i} = Stimulus{i}(:);
        Stimulus{i} = Stimulus{i}/max(abs(Stimulus{i}));
        StimDur(i) = length(Stimulus{i})/Fs;

        SortedStim = sort(Stimulus{i});
    else
        StimOnsets{i} = Info{2};
        StimOffsets{i} = Info{3};
        [TempData, Fs] = wavread(Info{1}{1});
        Silence = TempData(round(Fs*Info{2}):round(Fs*Info{3}));
        Silence = Silence(:);
    end
end

% Now put the different stimuli in a particular order. Only conditions that
% I follow is that total number of each stimulus should be equal to
% NumRepetitions

StimOrder = repmat([1:1:NumStimuli]', NumRepetitions, 1);
StimOrder = StimOrder(randperm(length(StimOrder)));
Intervals = (rand(size(StimOrder)) * (MaxInterval - MinInterval)) + MinInterval;

MaxNumStimPerFile = 6;

Fid = fopen([OutputFileName, '.Stimulus.List.txt'], 'w');
OutputStimulusFileIndex = 0;
OutputStimulus = [];
ElapsedTime = 0;
for i = 1:length(StimOrder),
    if (mod(i,MaxNumStimPerFile) == 1)
        OutputStimulusFileIndex = OutputStimulusFileIndex + 1;
        disp(['Output file #', num2str(OutputStimulusFileIndex), ' ...']);
        OutputStimFileName{OutputStimulusFileIndex} = [OutputFileName, '.', num2str(OutputStimulusFileIndex-1), '.wav'];
        if (~isempty(OutputStimulus))
            wavwrite(OutputStimulus(1:(length(OutputStimulus) - round(Fs*1)),:), Fs, 16, OutputStimFileName{OutputStimulusFileIndex});
        end
        OutputStimulus(1:(length(OutputStimulus) - round(Fs*1)),:) = [];
    end
    TempIntervalData = repmat(Silence(:), ceil(Intervals(i)) + 1, 1);
    TempIntervalData = TempIntervalData(1:round(Intervals(i)*Fs));
    if (length(Stimulus{StimOrder(i)}) < length(Silence))
        OutputStimulus = [OutputStimulus; [TempIntervalData TempIntervalData]; [OutputVolume*Stimulus{StimOrder(i)} Silence(1:length(Stimulus{StimOrder(i)}))]];
    else
        OtherSpeakerSilence = repmat(Silence(:), ceil(length(Stimulus{StimOrder(i)})/length(Silence)), 1);
        OutputStimulus = [OutputStimulus; [TempIntervalData TempIntervalData]; [OutputVolume*Stimulus{StimOrder(i)} OtherSpeakerSilence(1:length(Stimulus{StimOrder(i)}))]];
    end
    
    fprintf(Fid, 'Time %6.4fs: Stimulus - #%i - %s\n', ElapsedTime + Intervals(i), StimOrder(i), StimFileName{StimOrder(i)});
    ElapsedTime = ElapsedTime + Intervals(i) + StimDur(StimOrder(i));
    if (i == length(StimOrder))
        OutputStimFileName{OutputStimulusFileIndex} = [OutputFileName, '.', num2str(OutputStimulusFileIndex), '.wav'];
        OutputStimulus = [OutputStimulus; [Silence Silence]];
        if (~isempty(OutputStimulus))
            wavwrite(OutputStimulus(1:(length(OutputStimulus) - round(Fs*1)),:), Fs, 16, OutputStimFileName{OutputStimulusFileIndex});
        end
        OutputStimulus(1:(length(OutputStimulus) - round(Fs*1)),:) = [];
    end 
end

fclose(Fid);


disp('Finished');