function plotFF_combine(matfile1, matfile2, birdname_date);

% This script read two mat files for FF and plot and save them as
% an ascii file

% load matfiles
load (matfile1)
dir = [ffreq{:, 2}];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% added by Raghav
dir_duration = [ffreq{:,3}];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear ffreq

load (matfile2)
undir = [ffreq{:, 2}];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% added by Raghav
undir_duration = [ffreq{:,3}];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot FF
plot_notes4DS_UDS(' ', ' ', dir); hold on
plot_notes2(' ',' ', undir); axis auto;
title(['Fundamental Frequency', [birdname_date, ',  "', note, '"  (red: dir;  blue: undir)']])
% suptitle([birdname_date, ',  "', note, '"  (red: dir;  blue: undir)'])

% combine the data and save them
dir = dir'; undir = undir';
if length(dir)>length(undir)
    maxrow = length(dir)
elseif length(dir)<=length(undir)
    maxrow = length(undir);
end
ff_comb = zeros(maxrow, 2);
ff_comb(1:length(dir),1) = dir(:);
ff_comb(1:length(undir),2) = undir(:);

saved_file_name = ['FF_', note, '_comb']
save(saved_file_name, 'dir', 'undir')
save(saved_file_name, 'ff_comb', '-ascii', '-tabs')

fig_name =  ['FF_',note,'_comb.fig']
saveas(gcf,fig_name)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% added by Raghav
figure
plot(undir_duration,undir,'r.');
hold on;
plot(dir_duration,dir,'b.');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%