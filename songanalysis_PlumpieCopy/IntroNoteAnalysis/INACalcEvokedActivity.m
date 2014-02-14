function [MotifIntroNoteIndices, NonMotifIntroNoteIndices] = INACalcEvokedActivity(MotifIntroNoteIndices, NonMotifIntroNoteIndices, ASpikeData, FileLength)

Cols = ['rgbcmk'];
GaussianLen = 4;
Width = 0.005;
Fs = 10000;

XGauss = 1:1:(1 + round(2 * GaussianLen * Width * (Fs)));
XGauss = XGauss - (length(XGauss) + 1)/2;
GaussWin = (1/((Width * Fs) * sqrt(2 * pi))) * exp(-(XGauss.*XGauss)/(2 * (Width * Fs) * (Width * Fs)));

Time = (0:1/Fs:FileLength);
FR = zeros(size(Time));
Index = ceil((ASpikeData(1))*10000);
FR(1:Index-1) = 1/(ASpikeData(1));
for j = 1:length(ASpikeData),
    Index = ceil((ASpikeData(j))*10000);
    if (j == length(ASpikeData)),
        FR(Index:end) = 1/(Time(end) - ASpikeData(j));
    else
        NextIndex = ceil((ASpikeData(j+1))*10000);
        FR(Index:NextIndex-1) = 1/(ASpikeData(j+1) - ASpikeData(j));
    end
end

FR = conv(FR, GaussWin, 'same');

if (~isempty(MotifIntroNoteIndices))
    for i = 1:length(MotifIntroNoteIndices.NoofINs),
        if (MotifIntroNoteIndices.NoofINs{i} > 0)
            MotifIntroNoteIndices.SpontActivity{i}(:,1) = Time(find((Time >= (MotifIntroNoteIndices.Onsets{i}(1) - MotifIntroNoteIndices.VocalInterval{i})) & (Time <= MotifIntroNoteIndices.Onsets{i}(1))))';
            MotifIntroNoteIndices.SpontActivity{i}(:,2) = FR(find((Time >= (MotifIntroNoteIndices.Onsets{i}(1) - MotifIntroNoteIndices.VocalInterval{i})) & (Time <= MotifIntroNoteIndices.Onsets{i}(1))))';

            MotifIntroNoteIndices.SpontFiringRate{i} = length(find((ASpikeData >= (MotifIntroNoteIndices.Onsets{i}(1) - MotifIntroNoteIndices.VocalInterval{i})) & (ASpikeData <= MotifIntroNoteIndices.Onsets{i}(1))))/(MotifIntroNoteIndices.VocalInterval{i});
            MotifIntroNoteIndices.SpontSpikeTrain{i} = ASpikeData(find((ASpikeData >= (MotifIntroNoteIndices.Onsets{i}(1) - MotifIntroNoteIndices.VocalInterval{i})) & (ASpikeData <= MotifIntroNoteIndices.Onsets{i}(1))));
            MotifIntroNoteIndices.SpontNoofSpikes{i} = length(MotifIntroNoteIndices.SpontSpikeTrain{i});
        else
            MotifIntroNoteIndices.SpontActivity{i}(:,1) = Time(find((Time >= (MotifIntroNoteIndices.MotifOnsets{i}(1) - MotifIntroNoteIndices.VocalInterval{i})) & (Time <= MotifIntroNoteIndices.MotifOnsets{i}(1))))';
            MotifIntroNoteIndices.SpontActivity{i}(:,2) = FR(find((Time >= (MotifIntroNoteIndices.MotifOnsets{i}(1) - MotifIntroNoteIndices.VocalInterval{i})) & (Time <= MotifIntroNoteIndices.MotifOnsets{i}(1))))';

            MotifIntroNoteIndices.SpontFiringRate{i} = length(find((ASpikeData >= (MotifIntroNoteIndices.MotifOnsets{i}(1) - MotifIntroNoteIndices.VocalInterval{i})) & (ASpikeData <= MotifIntroNoteIndices.MotifOnsets{i}(1))))/(MotifIntroNoteIndices.VocalInterval{i});
            MotifIntroNoteIndices.SpontSpikeTrain{i} = ASpikeData(find((ASpikeData >= (MotifIntroNoteIndices.MotifOnsets{i}(1) - MotifIntroNoteIndices.VocalInterval{i})) & (ASpikeData <= MotifIntroNoteIndices.MotifOnsets{i}(1))));
            MotifIntroNoteIndices.SpontNoofSpikes{i} = length(MotifIntroNoteIndices.SpontSpikeTrain{i});
        end
    end
end
disp('Calculated evoked activity');