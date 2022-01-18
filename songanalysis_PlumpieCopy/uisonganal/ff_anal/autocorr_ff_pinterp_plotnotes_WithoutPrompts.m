function [FF] = autocorr_ff_pinterp_plotnotes_WithoutPrompts(Syllable, NoteSegment, Fs, FF_lo, FF_hi, StartPercent, EndPercent, PlotOption)


% Calculates FF of a note segment using autocorrelation and parabolic
% interpolation similar to Kao et al. 2005

Time = (1:1:length(Syllable))/Fs;
% plot the note
if (PlotOption == 1)
    figure;
    subplot(2,1,1);
    PlotSpectrogramInAxis_SongVar(Syllable, Time, Fs, subplot(2,1,1));

    % draw lines for the period of FF calculation
    FFStart = StartPercent * Time(end)/100;
    FFEnd = EndPercent * Time(end)/100;
    plot(ones(1,2)*FFStart, [300 8000], 'b-', 'LineWidth', 1.5);
    plot(ones(1,2)*FFEnd, [300 8000], 'b-', 'LineWidth', 1.5);
end

%calculate the auto-covariance of the song segment
%window in which it looks for the 1st peak after the 0th lag is 3-100 pts
[autocorr,lags] =xcov(NoteSegment);
%acwin=autocorr(length(note_segment)+3:length(note_segment)+70);  %for ~400-7000 Hz if Fs=32000 Hz
% if Fs=20000/ use range 3:40  
acwin=autocorr(length(NoteSegment)+round(Fs/FF_hi):length(NoteSegment)+round(Fs/FF_lo)); 

[max1,loc]=max(acwin);
if loc==1 | loc==length(acwin)
    peak=loc;   %parabolic interpolation req at least 3 points
else
%                 peak = loc;
    peak=pinterp([loc-1;loc;loc+1],[acwin(loc-1);acwin(loc);acwin(loc+1)]);
end     %peak values are real number
period=peak+(round(Fs/FF_hi) - 1);   %using window from 3rd-100th point after peak;
                %if peak is the 3rd point,peak=1; so period=peak+2                
FF=Fs/period;

% % SK modefied the follwoing part to limit highest ff value
% if FF > FF_hi
%    while 1
%        q = loc+1;
%        while q
%            if acwin(q-1)<acwin(q) & acwin(q) > acwin(q+1) 
%                loc3 = q;
%                break
%            end
%            q = q+1;
%        end
%        %[max1,loc2]=max(acwin(loc+5:length(acwin)));
%        %loc3 = loc2+loc+5;
%        peak = loc3;
%        peak=pinterp([loc3-1;loc3;loc3+1],[acwin(loc3-1);acwin(loc3);acwin(loc3+1)]);
%        period=peak+2;                               
%        FF=Fs/period;
%        if FF <= FF_hi
%            break
%        end
%        loc = loc3;
%    end
% end

if (PlotOption == 1)
    plot([FFStart FFEnd], FF*ones(1,2), 'g', 'LineWidth', 3)

    subplot(2,1,2);
    plot((1:1:length(acwin))+(round(Fs/FF_hi) - 1), acwin);
    hold on;
    plot(period, acwin(round(peak)), 'ks', 'MarkerFaceColor', 'k');
    subplot(2,1,1);
end