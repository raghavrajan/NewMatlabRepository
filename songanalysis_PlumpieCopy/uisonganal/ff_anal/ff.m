function [freq,clip] = ff(sample,Fs,iter)

% Estimate the fundamental frequency of a sample. 
% Current Methods
% 1 - One autocorr, find peak within a range and do parabolic interpolation
% 2 - Iterative autocorrelation until signal stabalizes, interpolate zeros and 
%     average to find period
% 3 - PSD autocorr - haven't really figured this one out yet
method = 1;

%calculate the auto-covariance of the song segment
%mimi's original with window 3-80 for ~400-7000 Hz
if method == 1
    autocorr=xcov(sample);
    acwin = autocorr(length(sample)+3:length(sample)+80);
    [y,loc]=max(acwin);
    %parabolic interpolation
    if  loc == 1 | loc == length(acwin)
        peak = loc; % can't do it without three points
    else
        peak = pinterp([loc-1;loc;loc+1],[acwin(loc-1);acwin(loc);acwin(loc+1)]);
    end
period=peak+2;
freq=Fs/period;
clip = autocorr(length(sample):end);

%iterative autocorrelation
elseif method == 2
    if nargin == 2
        iter = 120;
    end

    %performs iterative autocorrelation iter number of times
    xi = ixcov(sample,iter);

    %calculate average distance between zero crossings for first 150 points
    % (after that the quality of the numbers starts to degrade)
    xi=xi(1:150);
    zs = interpzeros(xi);
    d = diff(zs);
    if isempty(d) % no zero crossings
        period = Fs;
    else
        period = mean(d) * 2; % two zero crossings per period
    end
freq = Fs/period;
clip = xi;
% PSD autocorr - from mimi, unchanged
elseif method == 3
    %analyze the psd of selected song segment
    nfft=8192;
    window=8192;
    [Pxx,freq]=psd(sample,nfft,Fs,window);
    amp=sqrt(Pxx);
             
    %calculate the auto-covariance of the psd of the song segment
    autocorr=xcov(Pxx(1:1:500));
    [maxvalue,loc]=max(autocorr(530:750));
    loc=loc+29;
    freq=loc*(Fs/window);
end

