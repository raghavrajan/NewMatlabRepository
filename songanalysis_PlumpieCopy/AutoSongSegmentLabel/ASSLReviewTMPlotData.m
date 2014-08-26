function [SpecAxisLimits, LabelAxisLimits, AmpAxisLimits] = ASSLReviewTMPlotData(handles, Time, LogAmplitude)
% This function is called by all plots and this function then calls
% ASSLReviewPlotData

if (isfield(handles.ASSLReviewTMResults, 'SyllOnsets'))
    if (isfield(handles.ASSLReviewTMResults, 'SyllLabels'))
        if (get(handles.LoResHiResToggle, 'Value') == 1)
            [SpecAxisLimits, LabelAxisLimits, AmpAxisLimits] = ASSLReviewPlotData(handles.ASSLReviewTMResults.DirName, handles.ASSLReviewTMResults.FileName{handles.ASSLReviewTMResults.FileIndex}, handles.ASSLReviewTMResults.FileType, Time, LogAmplitude, handles.ReviewSpecAxis, handles.ReviewAmplitudeAxis, handles.ReviewLabelAxis, handles.ASSLReviewTMResults.Threshold{handles.ASSLReviewTMResults.FileIndex}, handles.ASSLReviewTMResults.SyllOnsets{handles.ASSLReviewTMResults.FileIndex}, handles.ASSLReviewTMResults.SyllOffsets{handles.ASSLReviewTMResults.FileIndex}, handles.ASSLReviewTMResults.SyllLabels{handles.ASSLReviewTMResults.FileIndex});
        else
            [SpecAxisLimits, LabelAxisLimits, AmpAxisLimits] = ASSLReviewPlotDataLowRes(handles.ASSLReviewTMResults.DirName, handles.ASSLReviewTMResults.FileName{handles.ASSLReviewTMResults.FileIndex}, handles.ASSLReviewTMResults.FileType, Time, LogAmplitude, handles.ReviewSpecAxis, handles.ReviewAmplitudeAxis, handles.ReviewLabelAxis, handles.ASSLReviewTMResults.Threshold{handles.ASSLReviewTMResults.FileIndex}, handles.ASSLReviewTMResults.SyllOnsets{handles.ASSLReviewTMResults.FileIndex}, handles.ASSLReviewTMResults.SyllOffsets{handles.ASSLReviewTMResults.FileIndex}, handles.ASSLReviewTMResults.SyllLabels{handles.ASSLReviewTMResults.FileIndex});
        end
    else
        if (handles.ASSLReviewTMResults.LoResHiRes == 1)
            [SpecAxisLimits, LabelAxisLimits, AmpAxisLimits] = ASSLReviewPlotData(handles.ASSLReviewTMResults.DirName, handles.ASSLReviewTMResults.FileName{handles.ASSLReviewTMResults.FileIndex}, handles.ASSLReviewTMResults.FileType, Time, LogAmplitude, handles.ReviewSpecAxis, handles.ReviewAmplitudeAxis, handles.ReviewLabelAxis, handles.ASSLReviewTMResults.Threshold{handles.ASSLReviewTMResults.FileIndex}, handles.ASSLReviewTMResults.SyllOnsets{handles.ASSLReviewTMResults.FileIndex}, handles.ASSLReviewTMResults.SyllOffsets{handles.ASSLReviewTMResults.FileIndex});
        else
            [handles.ASSLReviewTMResults.SpecAxisLimits, handles.ASSLReviewTMResults.LabelAxisLimits, handles.ASSLReviewTMResults.AmpAxisLimits] = ASSLReviewPlotDataLowRes(handles.ASSLReviewTMResults.DirName, handles.ASSLReviewTMResults.FileName{handles.ASSLReviewTMResults.FileIndex}, handles.ASSLReviewTMResults.FileType, Time, LogAmplitude, handles.ReviewSpecAxis, handles.ReviewAmplitudeAxis, handles.ReviewLabelAxis, handles.ASSLReviewTMResults.Threshold{handles.ASSLReviewTMResults.FileIndex}, handles.ASSLReviewTMResults.SyllOnsets{handles.ASSLReviewTMResults.FileIndex}, handles.ASSLReviewTMResults.SyllOffsets{handles.ASSLReviewTMResults.FileIndex});
        end
    end
else
    if (handles.ASSLReviewTMResults.LoResHiRes == 1)
        [handles.ASSLReviewTMResults.SpecAxisLimits, handles.ASSLReviewTMResults.LabelAxisLimits, handles.ASSLReviewTMResults.AmpAxisLimits] = ASSLReviewPlotData(handles.ASSLReviewTMResults.DirName, handles.ASSLReviewTMResults.FileName{handles.ASSLReviewTMResults.FileIndex}, handles.ASSLReviewTMResults.FileType, Time, LogAmplitude, handles.ReviewSpecAxis, handles.ReviewAmplitudeAxis, handles.ReviewLabelAxis);
    else
        [handles.ASSLReviewTMResults.SpecAxisLimits, handles.ASSLReviewTMResults.LabelAxisLimits, handles.ASSLReviewTMResults.AmpAxisLimits] = ASSLReviewPlotDataLowRes(handles.ASSLReviewTMResults.DirName, handles.ASSLReviewTMResults.FileName{handles.ASSLReviewTMResults.FileIndex}, handles.ASSLReviewTMResults.FileType, Time, LogAmplitude, handles.ReviewSpecAxis, handles.ReviewAmplitudeAxis, handles.ReviewLabelAxis);
    end
end