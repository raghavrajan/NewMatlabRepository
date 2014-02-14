function [LogAmplitude] = MA_CalculateLogAmplitudeAronovFee(RawData, Fs)

%==========================================================================
% Function for calculating the log amplitude of sound based on Aronov and
% Fee - J.Neurosci paper. Basically, song is filtered between 1000 and
% 4000Hz, squared and then it is smoothed with a 2.5ms window.
% This is then converted into log amplitude by doing 10*log10
%
% Raghav Rajan - 29th November 2013
%==========================================================================


F_High = 1000; % high pass filter cut-off
F_Low = 4000; % low pass filter cut-off

FilterForSong = fir1(80, [F_High*2/Fs F_Low*2/Fs], 'bandpass');
FiltSong = filtfilt(FilterForSong, 1, RawData);

SmoothWinSize = 2.5/1000;
 
Window = ones(round(SmoothWinSize*Fs), 1);
Window = Window/sum(Window);
LogAmplitude = 10*log10(conv(FiltSong.*FiltSong, Window, 'same'));

LogAmplitude = LogAmplitude(:)';