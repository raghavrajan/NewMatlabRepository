function [maxcorr_envelop] = xcorr_songenv(song1, song2, birdname);

% This script compares envelops of 2 songs using xcorr with time warping.
timewarp_low = 80; % lower tolerence for time warping (%)
timewarp_hi = 110; % higher tolerence for time warping (%)
Fs = 44100; % Sound frequency
F_low = 300; % Frequency for high pass fileter
F_hi = 8000; % Frequency for low pass fileter
filter_type = 'butter';
sm_win = 5; % smoothing window (msec)
increment = 5; % increment in xcorr analysis (msec)
dwnsample = 4.41; % downsampling rate of song (downsampled freq would be Fs/dwnsample)

% load songs
[rawsong1,Fs]=soundin('./', song1, 'w');
[rawsong2,Fs]=soundin('./', song2, 'w');

%filter songs
disp('filtering the songs...')
filtsong1=bandpass(rawsong1,Fs,F_low,F_hi,filter_type);
filtsong2=bandpass(rawsong2,Fs,F_low,F_hi,filter_type);

%calculate square of signal: proportional to power
squared_song1 = filtsong1.^2;
squared_song2 = filtsong2.^2;

% transform amplitude to dB scale
log_song1 = log10(squared_song1);
log_song2 = log10(squared_song2);

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
disp('downsampling the songs...')
for m = 1:fix(length(smsong1)/dwnsample)
    smsong1dwn(m) = mean(smsong1(round((m-1)*dwnsample+1):round(m*dwnsample)));
end
for m = 1:fix(length(smsong2)/dwnsample)
    smsong2dwn(m) = mean(smsong2(round((m-1)*dwnsample+1):round(m*dwnsample)));
end
smsong1 = smsong1dwn;
smsong2 = smsong2dwn;


  % measure cross-crrelation
  songlen1 = length(smsong1);
  songlen2 = length(smsong2);
  if songlen1>songlen2
      disp('ERROR!!  song2 must be longer than song1')
  end
  
timewarp = [timewarp_low:timewarp_hi];
increment = round(increment/1000*Fs);

for t=1:length(timewarp)
    song1X = [1:songlen1];
    song1X2 = song1X.*timewarp(t)./100;
    song1Xwarp = [2:fix(songlen1*timewarp(t)./100)];
    song1Ywarp = interp1(song1X2, smsong1, song1Xwarp);
    
    if length(song1Ywarp)>songlen2
        disp('too much timewarping!!  song1 > song2')
    elseif length(song1Ywarp)<=songlen2
        n=1;
        while n
            if (n-1)*increment+length(song1Ywarp) > length(smsong2)             
               break
           end
            smsong2win = smsong2((n-1)*increment+1:(n-1)*increment+length(song1Ywarp));
            coef(n) = max(xcov(song1Ywarp, smsong2win, round(increment/2), 'coeff')); % calculate corr.coef after subtracting mean of each song
            n = n+1;
        end  
         xcorr_envelop(t) = max(coef);
         display = ['coefficient = ',num2str(xcorr_envelop(t)),' (',num2str(timewarp(t)),'% time-warped song)'];
         disp(display)
         clear coef
     end
 end
 
 maxcorr_envelop = max(xcorr_envelop)
 
 %save the results
 saved_file_name = ['XcorrEnv_', song1, '_', song2]
 save(saved_file_name,'maxcorr_envelop','xcorr_envelop','timewarp_low', 'timewarp_hi')
 
 % plot the songs and coeficients
 figure;
 subplot(3,1,1); plot(smsong1); title('reference song');axis([0 length(smsong2) max(smsong1)*(-0.1) max(smsong1)*1.1])
 subplot(3,1,2); plot(smsong2); title('test song');axis([0 length(smsong2) max(smsong2)*(-0.1) max(smsong2)*1.1])
 subplot(3,1,3); plot(timewarp, xcorr_envelop); title('correlation coeficient with different time warping of reference song (% change)');
 axis([min(timewarp) max(timewarp) 0 1]); 
 xlabel('% change from original song1'); ylabel('corr. coefficient')
 ttl = ['cross-correlation between ',song1,' & ',song2,' (',birdname,')' ];
 suptitle(ttl)
 
 fig_name = [saved_file_name,'.fig']
saveas(gcf,fig_name)
close all