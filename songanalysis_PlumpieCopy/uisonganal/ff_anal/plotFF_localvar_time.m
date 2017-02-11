function plotFF_localvar_time(note);

% This script reads FF_note_localvar3.mat file and calculates SD of %
% difference using sliding time windows to see how local variability changes over time.
% Also, calculates a correlation coefficient between SD time couces of dir
% and that of undir.

matfile = ['FF_',note,'_localvar3.mat']
load(matfile);

% calculate local variance using a sliding time window

windur = 30 % window duration (min)
increment = 1 % window increment (min)
windur = windur/60;
increment = increment/60;

windstart_dir = fix(min(dirX)-windur);
n = 1;
i = windstart_dir;
while i
    localdata_dir = find(dirX >= i & dirX <= i+windur);
    if length(localdata_dir)>4
        localdata_dirY = deltaF_dir(localdata_dir);
        localstd_dir(n) = std(localdata_dirY);
        localdata_dirX(n) = i + windur/2;
        n = n+1;
    end
    i = i+increment;
    if i > max(dirX)
        break
    end
end

windstart_undir = fix(min(undirX)-windur);
n = 1;
i = windstart_undir;
while i
    localdata_undir = find(undirX >= i & undirX <= i+windur);
    if length(localdata_undir)>4
        localdata_undirY = deltaF_undir(localdata_undir);
        localstd_undir(n) = std(localdata_undirY);
        localdata_undirX(n) = i + windur/2;
        n = n+1;
    end
    i = i+increment;
    if i > max(undirX)
        break
    end
end

% plot
figure
subplot(2,1,1)
plot(dirX,deltaF_dir,'or'); hold on
plot(undirX,deltaF_undir,'ob')
title('% difference from local mean (calculated 30 min window around each data point)')
xlabel('time'); ylabel('% difference')
zeroX = [fix(min([min(dirX), min(undirX)])) ceil(max([max(dirX), max(undirX)]))];
zeroY = [0 0];
plot(zeroX, zeroY, '--k')

subplot(2,1,2)
plot(localdata_dirX,localstd_dir,'or'); hold on
plot(localdata_undirX,localstd_undir,'ob'); hold on
title('SD of % difference (calculated 30 min sliding windows with 1 min increments)')

FF_time_matfile = ['FF_',note,'_time.mat'];
suptitle([FF_time_matfile,  ' (red: dir;  blue: undir)'])

% Calculate correlation coefficient of SDs between dir & undir
coefXstart = max(min(localdata_dirX), min(localdata_undirX)); 
coefXend = min(max(localdata_dirX), max(localdata_undirX));% Find overlapping parts between dir & undir
coefX = [coefXstart:increment:coefXend];

coefstd_dir = interp1(localdata_dirX,localstd_dir,coefX); % interporate the data points
coefstd_undir = interp1(localdata_undirX,localstd_undir,coefX);

plot(coefX,coefstd_dir,'+k'); hold on
plot(coefX,coefstd_undir,'+k'); hold on

corr_coef = corrcoef(coefstd_dir, coefstd_undir);

corr_coef_text = ['corr.coef = ',num2str(corr_coef(1,2))]
x =     min([min(localdata_dirX), min(localdata_undirX)]);
text(x, max(coefstd_undir)*1.1, corr_coef_text)

% save the variables
saved_file_name = ['FF_', note, '_localvar_time']
save(saved_file_name,'localstd_dir','localdata_dirX','localstd_undir','localdata_undirX',...
    'coefstd_dir', 'coefstd_undir','coefX','corr_coef')

% save the figure
fig_name =  ['FF_',note,'_localvar_time.fig']
saveas(gcf,fig_name)