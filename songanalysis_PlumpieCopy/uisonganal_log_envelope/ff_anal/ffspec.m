function ffspec(song, fs)

%plots the spectrogram of the song
%user clicks what he/she thinks is the ff; will plot the autocorrelation of the
%signal and tell the ff estimate, the user's guess, and the ratio between the two
win = 250;

subplot(2,1,1)
specgram(song,300,fs,280,250);
%axis([0 length(song)/fs 0 8000])
set(gca,'ylim',[0 8000])
caxis([40 100])
colormap(jet)

while 1
    [x,y] = ginput(1);
    
    pos = fix(x * fs);
    sample = song(pos-win:pos+win);
    subplot(2,1,2)
    [freq,clip] = ff(sample,fs);
    plot(clip)
    set(gca, 'xlim', [0 150]);
    freqText = xlabel(['ff estimate: ' num2str(freq) '   y: ' num2str(y) '   ff/y: ' num2str(freq/y)]);
end