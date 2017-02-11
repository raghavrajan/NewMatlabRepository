function ffplot(song,fs)

%plots the spectrogram and 1000 random points of the frequency estimate (using ff.m)

specgram(song,300,fs,280,250);
caxis([-30 10])
colormap(flipud(gray))
axis([0 length(song)/fs 0 10000])

hold on
pause(.2)
for i = 1:1000
    i = ceil(rand*(length(song)-600));
    ffs = ff(song(i:i+500),fs);
    p=plot(i/fs,ffs,'r.');
    set(p,'markersize',6) 
end

