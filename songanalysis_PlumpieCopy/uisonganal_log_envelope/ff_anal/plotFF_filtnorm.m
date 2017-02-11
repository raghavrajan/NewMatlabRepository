function plotFF_filtnorm (FF_time_matfile, note);



% frequency of low pass filter
F_low = 0.001
Fs = 10000

% load input files
load (FF_time_matfile)

% sort the FF values for dir and undir
m =1; p =1;
for k=1:length(allFF)
    if length(strmatch('dir', allFF(k,1)))>0
        dirX(m) = [allFF{k,3}];
        dirY(m) = [allFF{k,2}];
        m=m+1;
    end
    if length(strmatch('undir', allFF(k,1)))>0
        undirX(p) = [allFF{k,3}];
        undirY(p) = [allFF{k,2}];
        p=p+1;
    end
end

% modify X for interpolation
dirX2= dirX;
for x = 1:length(dirX)-1
    if dirX(x+1) == dirX(x)
        dirX2(x+1) = dirX2(x)+0.0001;
    end
end
undirX2= undirX;
for x = 1:length(undirX)-1
    if undirX(x+1) == undirX(x)
        undirX2(x+1) = undirX2(x)+0.0001;
    end
end

% interpolation
dirX_int = [min(dirX2):0.0001:max(dirX2)];
dirY_int = interp1(dirX2, dirY,dirX_int);
undirX_int = [min(undirX2):0.0001:max(undirX2)];
undirY_int = interp1(undirX2, undirY ,undirX_int);

% filter the signals
nfilt = 20000
lowpass = fir1(nfilt, F_low*2/Fs, 'low');
dirY_filt = filtfilt(lowpass, 1, dirY_int); 







% fit with a regression curve
dirFit = polyfit(dirX,dirY,n);
undirFit = polyfit(undirX,undirY,n);

dirFitY = polyval(dirFit,dirX);
undirFitY = polyval(undirFit,undirX);

% plot FF varlues and regression curves
figure
subplot(2,1,1)
plot(dirX,dirY,'or'); hold on
plot(dirX,dirFitY,'--r'); hold on
plot(undirX,undirY,'ob'); hold on
plot(undirX,undirFitY,'--b')
title(['Regression curves (degree: ', int2str(n), '),  ', FF_time_matfile,  ',  Red: dir,  Blue: undir'])
xlabel('time'); ylabel('FF (Hz)')

% plot difference between FF value and regression curve
dirYdif = dirY-dirFitY;
undirYdif = undirY-undirFitY;
subplot(2,1,2)
plot(dirX,dirYdif,'or'); hold on
plot(undirX,undirYdif,'ob'); hold on
axis auto
xall = [allFF{:,3}]; y = xall*0;
plot(xall,y,'--k')
xlabel('time'); ylabel('difference from regression line (Hz)')

%plot the normalized values
plot_notes4DS_UDS(' ', ' ', dirYdif); hold on
plot_notes2(' ',' ', undirYdif); axis auto;
title(['Difference from regression curve (degree: ', int2str(n), '),  ', FF_time_matfile,'  Red: dir,  Blue: undir'])

% save the FF difference
saved_file_name = ['FF_', note, '_regression2']
save(saved_file_name, 'dirYdif', 'undirYdif')