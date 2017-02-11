function [MotifIntroNoteIndices, NonMotifIntroNoteIndices] = INACalcEvokedActivity2(MotifIntroNoteIndices, NonMotifIntroNoteIndices, ASpikeData, Labels, FileLength)

Cols = ['rgbcmk'];
GaussianLen = 4;
Width = 0.005;
Fs = 10000;
PreDuration = 0.06;

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
            for j = 1:(MotifIntroNoteIndices.NoofINs{i}),
                Indices = (find((ASpikeData > (MotifIntroNoteIndices.Onsets{i}(j) - PreDuration)) & (ASpikeData <= MotifIntroNoteIndices.Onsets{i}(j))));
                if (~isempty(Indices))
                    MotifIntroNoteIndices.PreSpikeTrain{i}{j} = ASpikeData(Indices) - MotifIntroNoteIndices.Onsets{i}(j);
                    MotifIntroNoteIndices.PreFiringRate{i}(j,1) = length(ASpikeData(Indices))/PreDuration;
                    MotifIntroNoteIndices.PreSmoothedFiringRate{i}{j} = FR(find((Time > (MotifIntroNoteIndices.Onsets{i}(j) - PreDuration)) & (Time <= MotifIntroNoteIndices.Onsets{i}(j))))';
                else
                    MotifIntroNoteIndices.PreFiringRate{i}(j,1) = 0;
                end
                
                Indices = (find((ASpikeData > (MotifIntroNoteIndices.Onsets{i}(j))) & (ASpikeData <= MotifIntroNoteIndices.Offsets{i}(j))));
                if (~isempty(Indices))
                    MotifIntroNoteIndices.DurSpikeTrain{i}{j} = ASpikeData(Indices) - MotifIntroNoteIndices.Onsets{i}(j);
                    MotifIntroNoteIndices.DurFiringRate{i}(j,1) = length(ASpikeData(Indices))/(MotifIntroNoteIndices.Durs{i}(j));
                    MotifIntroNoteIndices.DurSmoothedFiringRate{i}{j} = FR(find((Time > (MotifIntroNoteIndices.Onsets{i}(j))) & (Time <= MotifIntroNoteIndices.Offsets{i}(j))))';
                else
                    MotifIntroNoteIndices.DurFiringRate{i}(j,1) = 0;
                end
                MotifIntroNoteIndices.NextSyll{i}(j,1) = Labels(MotifIntroNoteIndices.Indices{i}(j) + 1);
            end
        end
    end
end
disp('Calculated evoked activity');