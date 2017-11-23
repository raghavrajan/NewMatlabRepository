function [] = CompareHFMic_MixerMic_FF(HFMicDir, HFMicFileList, MixerMicDir, MixerMicFileList, FileType, SyllLabel, FF_Low, FF_High, PlotOption, FFMeasureDur, InputPercentFromStart)

PresentDir = pwd;

cd(HFMicDir);
[HFMicffreq] = autocorr_ff_pinterp_plotnotes_yincomparisonWithoutPrompts(HFMicDir, HFMicFileList, SyllLabel, FF_Low, FF_High, FileType, PlotOption, FFMeasureDur, InputPercentFromStart);

cd(MixerMicDir);
[MixerMicffreq] = autocorr_ff_pinterp_plotnotes_yincomparisonWithoutPrompts(MixerMicDir, MixerMicFileList, SyllLabel, FF_Low, FF_High, FileType, PlotOption, FFMeasureDur, InputPercentFromStart);

disp('Finished');