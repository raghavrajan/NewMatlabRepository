function plotFF_localvar2 (note);

% This script calculate CV of FF values using a time window around each data point.

windur = 30 % window duration (min)

% load input files
FF_time_matfile = ['FF_',note,'_time.mat']
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

% calculate local variance around each data point

windur = windur/60;
n=1;
for i = 1:length(dirX)
    localdata_dir = find(dirX >= dirX(i)-(windur/2) & dirX <= dirX(i)+(windur/2));
    if length(localdata_dir)>4
        localdata_dirY = dirY(localdata_dir);
        localstd_dir(n) = std(localdata_dirY);
        localmean_dir(n) = mean(localdata_dirY);
        localcv_dir(n) = localstd_dir(n)/localmean_dir(n);
        localdata_dirX(n) = dirX(i);
        n=n+1;
    end
end
n=1;
for i = 1:length(undirX)
    localdata_undir = find(undirX >= undirX(i)-(windur/2) & undirX <= undirX(i)+(windur/2));
    if length(localdata_undir)>4
        localdata_undirY = undirY(localdata_undir);
        localstd_undir(n) = std(localdata_undirY);
        localmean_undir(n) = mean(localdata_undirY);
        localcv_undir(n) = localstd_undir(n)/localmean_undir(n);
        localdata_undirX(n) = undirX(i);
        n=n+1;
    end
end


% plot FF varlues and CV
figure
subplot(2,1,1)
plot(dirX,dirY,'or'); hold on
plot(undirX,undirY,'ob'); hold on
title('Fundamental Frequency')
ylabel('FF (Hz)')

subplot(2,1,2)
plot(localdata_dirX,localcv_dir,'or'); hold on
plot(localdata_undirX,localcv_undir,'ob')
title('C.V. of FF (30 min window around each data point)')
xlabel('time'); ylabel('C.V')

suptitle([FF_time_matfile,  ' (red: dir;  blue: undir)'])

% save the FF difference
saved_file_name = ['FF_', note, '_localvar2']
save(saved_file_name,'localstd_dir','localmean_dir','localcv_dir','localdata_dirX',...
    'localstd_undir','localmean_undir','localcv_undir','localdata_undirX')
clear