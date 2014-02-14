function [maxcorr_envelop] = xcorr_songenv(song1, song2, birdname,template_file,template_no);

% This script compares envelops of 2 songs using xcorr with time warping.


%timewarp_low = input ('How much do you want to compress the songs (% of original song)?  ');
%timewarp_hi = input ('How much do you want to stretch the songs (% of original song)?  ');
timewarp_low = 80; % lower tolerence for time warping (%)
timewarp_hi = 120; % higher tolerence for time warping (%)

%%NOTE:  sampling rate is hard coded
Fs = 44100; % Sound frequency
%Fs=32000;

%filter the song between 300Hz and 8 kHz
F_low = 300; % Frequency for high pass filter
F_hi = 8000; % Frequency for low pass filter

%bandpass uses the butterworth filter
%filter_type = 'butter';

sm_win = 5; % smoothing window (msec)
increment = 5; % increment in xcorr analysis (msec)
%dwnsample = 3.2;
dwnsample = 4.41; % downsampling rate of song (downsampled freq would be Fs/dwnsample;
                  %in this case = 10K)

% load songs; assumes these are .wav files
%[rawsong1,Fs]=soundin('./', song1, 'w');
%song1 is the template; pulled out from the raw data
rawsong1=song1;    %motif template
%[rawsong2,Fs]=soundin('./', song2, 'obs0r');   %load the bout if observer
[rawsong2,Fs]=soundin('./', song2, 'w');   %load the bout



%filter songs
disp('filtering the songs...')
%filtsong1=bandpass(rawsong1,Fs,F_low,F_hi,filter_type);
%filtsong2=bandpass(rawsong2,Fs,F_low,F_hi,filter_type);
%modified in july 2013 b/c bandpass.m uses a butterworth filter; don't need to specify the filter type
filtsong1=bandpass(rawsong1,Fs,F_low,F_hi);
filtsong2=bandpass(rawsong2,Fs,F_low,F_hi);

%calculate square of signal: proportional to power
squared_song1 = filtsong1.^2;
squared_song2 = filtsong2.^2;

% transform amplitude to dB scale;  modified sk's code to add the factor of 10
log_song1 = 10*log10(squared_song1);
log_song2 = 10*log10(squared_song2);

%%july 2013 - mk
%%compare this with matlab function mag2db; gives similar values
%note:  from gk, should use 20log10(sigma) for db; or can use 10log10(x^2)
%mag2dbm.m
%y(y<0) = NaN;
%ydb = 20*log10(y);
%filtsong1_db=mag2db(filtsong1);
%filtsong2_db=mag2db(filtsong2);


%smooth the rectified song
disp('smoothing the songs...')
len=round(Fs*sm_win/1000);                      
h=ones(1,len)/len;
  
smooth1=conv(h, log_song1);
offset1=round((length(smooth1)-length(filtsong1))/2); %get rid of convolution induced offset
smsong1=smooth1(1+offset1:length(filtsong1)+offset1);
  
smooth2=conv(h, log_song2);
offset2=round((length(smooth2)-length(filtsong2))/2); %get rid of convolution induced offset
smsong2=smooth2(1+offset2:length(filtsong2)+offset2);

% downsample the songs to calculate faster
% so, here, we calculate the mean of x no. of points when downsampling; does not use resample.m
disp('downsampling the songs...')
%fix rounds to nearest integer towards 0
for m = 1:fix(length(smsong1)/dwnsample)
    smsong1dwn(m) = mean(smsong1(round((m-1)*dwnsample+1):round(m*dwnsample)));
end
for m = 1:fix(length(smsong2)/dwnsample)
    smsong2dwn(m) = mean(smsong2(round((m-1)*dwnsample+1):round(m*dwnsample)));
end
%smsong1 = smsong1dwn;
%smsong2 = smsong2dwn;

% measure cross-correlation of the smoothed, downsampled signal
  songlen1 = length(smsong1dwn)
  songlen2 = length(smsong2dwn)
  if songlen1>songlen2
      disp('ERROR!!  song2 must be longer than song1')
  end
  
timewarp = [timewarp_low:timewarp_hi];
increment = round(increment/1000*Fs);  %convert from ms to samples

for t=1:length(timewarp)
    song1X = [1:songlen1];
    song1X2 = song1X.*timewarp(t)./100;
    song1Xwarp = [2:fix(songlen1*timewarp(t)./100)];
    song1Ywarp = interp1(song1X2, smsong1dwn, song1Xwarp);
    
    if length(song1Ywarp)>songlen2
        disp('too much timewarping!!  song1 > song2')
    elseif length(song1Ywarp)<=songlen2
        n=1;
        while n
            if (n-1)*increment+length(song1Ywarp) > length(smsong2dwn)             
               break
           end
           
            smsong2win = smsong2dwn((n-1)*increment+1:(n-1)*increment+length(song1Ywarp));  %move along the song with a window equal to length of template
            coef(n) = max(xcov(song1Ywarp, smsong2win, round(increment/2), 'coeff')); % calculate corr.coef after subtracting mean of each song (xcorr does not remove the mean)
                                                                                      % song1Ywarp is the warped template; smsong2win is the window in song2
            n = n+1;
        end  
         %  so, coef has the max xcov for each window
         xcorr_envelop(t) = max(coef);
         display = ['coefficient = ',num2str(xcorr_envelop(t)),' (',num2str(timewarp(t)),'% time-warped song)'];
         disp(display)
         clear coef
     end
 end
 
maxcorr_envelop = max(xcorr_envelop)
maxcorr_timewarp_idx=find(xcorr_envelop==max(xcorr_envelop));
maxcorr_timewarp=timewarp(maxcorr_timewarp_idx);
display=['max coefficient with timewarp = ',(num2str(maxcorr_timewarp)),'% of the template'];
disp(display)
 
 
 %save the results
 saved_file_name = ['XcorrEnv_template', num2str(template_no),'_', song2,'.mat']
 save(saved_file_name,'maxcorr_envelop','xcorr_envelop','timewarp_low', 'timewarp_hi','maxcorr_timewarp')
 
 % plot the songs and coeficients
 figure;
 subplot(3,1,1); plot(smsong1dwn); title('reference song');axis([0 length(smsong2dwn) max(smsong1dwn)*(-0.1) max(smsong1dwn)*1.1])
 subplot(3,1,2); plot(smsong2dwn); title('test song');axis([0 length(smsong2dwn) max(smsong2dwn)*(-0.1) max(smsong2dwn)*1.1])
 subplot(3,1,3); plot(timewarp, xcorr_envelop); 
 l=line([maxcorr_timewarp maxcorr_timewarp ],[0 1]);
 set(l,'color','r')
 hold on
 t=title(['correlation coeficient with time warping of reference song (% change); max at ',num2str(maxcorr_timewarp),'%']);
 set(t,'fontsize',8)
 orient tall
 axis([min(timewarp) max(timewarp) 0 1]); 
 xlabel('% change from original song'); ylabel('corr. coefficient')
 ttl = ['xcorrelation betw ',template_file,' & ',song2,' (',birdname,')' ];
 suptitle(ttl)

fig_name = ['XcorrEnv_template', num2str(template_no),'_', song2,'.fig']
saveas(gcf,fig_name)
 
clear song1 song2 song1Ywarp smsong2win xcorr_envelop

close all