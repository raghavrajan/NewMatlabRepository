function [] = ASSLPlotSyllAmplitudeWaveforms(SyllData, Condition)

switch Condition
    case 'NormalOnset'
        for i = 1:length(SyllData.Amplitude),
            plot(SyllData.Time{i}, SyllData.Amplitude{i}, 'b');
            hold on;
        end
        
    case 'AdjustedOnset'
        for i = 1:length(SyllData.Amplitude),
            plot(SyllData.Time{i} - SyllData.Adjust(i,1), SyllData.Amplitude{i}, 'b');
            hold on;
        end
        
   case 'NormalOffset'
        for i = 1:length(SyllData.Amplitude),
            plot(SyllData.OffsetAlignedTime{i}, SyllData.Amplitude{i}, 'b');
            hold on;
        end
        
    case 'AdjustedOffset'
        for i = 1:length(SyllData.Amplitude),
            plot(SyllData.OffsetAlignedTime{i} - SyllData.Adjust(i,2), SyllData.Amplitude{i}, 'b');
            hold on;
        end
end
       
axis tight;
temp = axis;
plot([0 0], [temp(3) temp(4)], 'k--');

