function [SpecAxisLimits, LabelAxisLimits, AmpAxisLimits] = ASSLReviewTMPlotData(handles, Time, LogAmplitude)
% This function is called by all plots and this function then calls
% ASSLReviewPlotData


FieldNames = fieldnames(handles);
Index = find(cellfun(@length, strfind(FieldNames, 'ASSL')));

DataStruct = eval(['handles.', FieldNames{Index}]);

if (isfield(DataStruct, 'SyllOnsets'))
    if (isfield(DataStruct, 'SyllLabels'))
        if (get(handles.LoResHiResToggle, 'Value') == 1)
            [SpecAxisLimits, LabelAxisLimits, AmpAxisLimits] = ASSLReviewPlotData(DataStruct.FileDir{DataStruct.FileIndex}, DataStruct.FileName{DataStruct.FileIndex}, DataStruct.FileType, Time, LogAmplitude, handles.ReviewSpecAxis, handles.ReviewAmplitudeAxis, handles.ReviewLabelAxis, DataStruct.Threshold{DataStruct.FileIndex}, DataStruct.SyllOnsets{DataStruct.FileIndex}, DataStruct.SyllOffsets{DataStruct.FileIndex}, DataStruct.SyllLabels{DataStruct.FileIndex});
        else
            [SpecAxisLimits, LabelAxisLimits, AmpAxisLimits] = ASSLReviewPlotDataLowRes(DataStruct.FileDir{DataStruct.FileIndex}, DataStruct.FileName{DataStruct.FileIndex}, DataStruct.FileType, Time, LogAmplitude, handles.ReviewSpecAxis, handles.ReviewAmplitudeAxis, handles.ReviewLabelAxis, DataStruct.Threshold{DataStruct.FileIndex}, DataStruct.SyllOnsets{DataStruct.FileIndex}, DataStruct.SyllOffsets{DataStruct.FileIndex}, DataStruct.SyllLabels{DataStruct.FileIndex});
        end
    else
        if (DataStruct.LoResHiRes == 1)
            [SpecAxisLimits, LabelAxisLimits, AmpAxisLimits] = ASSLReviewPlotData(DataStruct.FileDir{DataStruct.FileIndex}, DataStruct.FileName{DataStruct.FileIndex}, DataStruct.FileType, Time, LogAmplitude, handles.ReviewSpecAxis, handles.ReviewAmplitudeAxis, handles.ReviewLabelAxis, DataStruct.Threshold{DataStruct.FileIndex}, DataStruct.SyllOnsets{DataStruct.FileIndex}, DataStruct.SyllOffsets{DataStruct.FileIndex});
        else
            [SpecAxisLimits, LabelAxisLimits, AmpAxisLimits] = ASSLReviewPlotDataLowRes(DataStruct.FileDir{DataStruct.FileIndex}, DataStruct.FileName{DataStruct.FileIndex}, DataStruct.FileType, Time, LogAmplitude, handles.ReviewSpecAxis, handles.ReviewAmplitudeAxis, handles.ReviewLabelAxis, DataStruct.Threshold{DataStruct.FileIndex}, DataStruct.SyllOnsets{DataStruct.FileIndex}, DataStruct.SyllOffsets{DataStruct.FileIndex});
        end
    end
else
    if (DataStruct.LoResHiRes == 1)
        [SpecAxisLimits, LabelAxisLimits, AmpAxisLimits] = ASSLReviewPlotData(DataStruct.FileDir{DataStruct.FileIndex}, DataStruct.FileName{DataStruct.FileIndex}, DataStruct.FileType, Time, LogAmplitude, handles.ReviewSpecAxis, handles.ReviewAmplitudeAxis, handles.ReviewLabelAxis);
    else
        [SpecAxisLimits, LabelAxisLimits, AmpAxisLimits] = ASSLReviewPlotDataLowRes(DataStruct.FileDir{DataStruct.FileIndex}, DataStruct.FileName{DataStruct.FileIndex}, DataStruct.FileType, Time, LogAmplitude, handles.ReviewSpecAxis, handles.ReviewAmplitudeAxis, handles.ReviewLabelAxis);
    end
end

if (handles.ASSLReviewTMResults.ScaleBarToggleValue == 1)
    if (isfield(handles.ASSLReviewTMResults, 'SyllOnsets'))
        Onsets = [handles.ASSLReviewTMResults.SyllOnsets{handles.ASSLReviewTMResults.FileIndex}(:); handles.ASSLReviewTMResults.FileDur{handles.ASSLReviewTMResults.FileIndex}*1000];
        Offsets = [0; handles.ASSLReviewTMResults.SyllOffsets{handles.ASSLReviewTMResults.FileIndex}(:)];
        Intervals = Onsets - Offsets;
        LongIntervals = find(Intervals >= handles.ASSLReviewTMResults.LongIntervalCutoff*1000);
        if (~isempty(LongIntervals))
            for i = LongIntervals(:)',
                axes(handles.ReviewSpecAxis);
                hold on;
                PatchXValues = [Offsets(i) Offsets(i) Onsets(i) Onsets(i)]/1000;
                patch(PatchXValues, [handles.ASSLReviewTMResults.SpecAxisLimits(3:4) fliplr(handles.ASSLReviewTMResults.SpecAxisLimits(3:4))], 'g', 'FaceColor', 'g', 'EdgeColor', 'none', 'FaceAlpha', 0.15);

                axes(handles.ReviewAmplitudeAxis);
                hold on;
                patch(PatchXValues, [handles.ASSLReviewTMResults.AmpAxisLimits(3:4) fliplr(handles.ASSLReviewTMResults.AmpAxisLimits(3:4))], 'g', 'FaceColor', 'g', 'EdgeColor', 'none', 'FaceAlpha', 0.15);
            end
        end
    end
end
ASSLAdjustSpectrogram(handles.ASSLReviewTMResults.SpectFloor, handles.ASSLReviewTMResults.SpectCeil, handles.ASSLReviewTMResults.Brightness, 'hot', 'classic', handles.ReviewSpecAxis);