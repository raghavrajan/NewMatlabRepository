function plotFF_summary(directory_pre1, directory_pre2, directory_post1d, directory_post2wk, directory_post8wk, bird, note);

% This script read FF file (e.g. FF_a_comb) and plot the mean, SD, and CV.
% Also, read normalized FF file (e.g. FF_a_meansubtract) and plot SD.


%%%%%% read FF files and calculate mean etc. %%%%%%

pre1 = [directory_pre1,'/FF_',note,'_comb.mat'];
pre2 = [directory_pre2,'/FF_',note,'_comb.mat'];
post1d = [directory_post1d,'/FF_',note,'_comb.mat'];
post2wk = [directory_post2wk,'/FF_',note,'_comb.mat'];
post8wk = [directory_post8wk,'/FF_',note,'_comb.mat'];

% load files for FF
temp = load (pre1);
dir_pre1 = getfield(temp,'dir'); undir_pre1 = getfield(temp,'undir');
temp = load (pre2);
dir_pre2 = getfield(temp,'dir'); undir_pre2 = getfield(temp,'undir');
temp = load (post1d);
dir_post1d = getfield(temp,'dir'); undir_post1d = getfield(temp,'undir');
temp = load (post2wk);
dir_post2wk = getfield(temp,'dir'); undir_post2wk = getfield(temp,'undir');
%temp = load (post4wk);
%dir_post4wk = getfield(temp,'dir'); undir_post4wk = getfield(temp,'undir');
temp = load (post8wk);
dir_post8wk = getfield(temp,'dir'); undir_post8wk = getfield(temp,'undir');

% calculate mean for FF
dir_pre1_mean = mean(dir_pre1);
undir_pre1_mean = mean(undir_pre1);

dir_pre2_mean = mean(dir_pre2);
undir_pre2_mean = mean(undir_pre2);

dir_post1d_mean = mean(dir_post1d);
undir_post1d_mean = mean(undir_post1d);

dir_post2wk_mean = mean(dir_post2wk);
undir_post2wk_mean = mean(undir_post2wk);

%dir_post4wk_mean = mean(dir_post4wk);
%undir_post4wk_mean = mean(undir_post4wk);

dir_post8wk_mean = mean(dir_post8wk);
undir_post8wk_mean = mean(undir_post8wk);

% calculate SD for FF
dir_pre1_std = std(dir_pre1);
undir_pre1_std = std(undir_pre1);

dir_pre2_std = std(dir_pre2);
undir_pre2_std = std(undir_pre2);

dir_post1d_std = std(dir_post1d);
undir_post1d_std = std(undir_post1d);

dir_post2wk_std = std(dir_post2wk);
undir_post2wk_std = std(undir_post2wk);

%dir_post4wk_std = std(dir_post4wk);
%undir_post4wk_std = std(undir_post4wk);

dir_post8wk_std = std(dir_post8wk);
undir_post8wk_std = std(undir_post8wk);

% calculate SEM for FF
dir_pre1_sem = std(dir_pre1)/(length(dir_pre1))^(1/2);
undir_pre1_sem = std(undir_pre1)/(length(undir_pre1))^(1/2);

dir_pre2_sem = std(dir_pre2)/(length(dir_pre2))^(1/2);
undir_pre2_sem = std(undir_pre2)/(length(undir_pre2))^(1/2);

dir_post1d_sem = std(dir_post1d)/(length(dir_post1d))^(1/2);
undir_post1d_sem = std(undir_post1d)/(length(undir_post1d))^(1/2);

dir_post2wk_sem = std(dir_post2wk)/(length(dir_post2wk))^(1/2);
undir_post2wk_sem = std(undir_post2wk)/(length(undir_post2wk))^(1/2);

%dir_post4wk_sem = std(dir_post4wk)/(length(dir_post4wk))^(1/2);
%undir_post4wk_sem = std(undir_post4wk)/(length(undir_post4wk))^(1/2);

dir_post8wk_sem = std(dir_post8wk)/(length(dir_post8wk))^(1/2);
undir_post8wk_sem = std(undir_post8wk)/(length(undir_post8wk))^(1/2);

% calculate CV for FF
dir_pre1_cv = 100*dir_pre1_std/dir_pre1_mean;
undir_pre1_cv = 100*undir_pre1_std/undir_pre1_mean;

dir_pre2_cv = 100*dir_pre2_std/dir_pre2_mean;
undir_pre2_cv = 100*undir_pre2_std/undir_pre2_mean;

dir_post1d_cv = 100*dir_post1d_std/dir_post1d_mean;
undir_post1d_cv = 100*undir_post1d_std/undir_post1d_mean;

dir_post2wk_cv = 100*dir_post2wk_std/dir_post2wk_mean;
undir_post2wk_cv = 100*undir_post2wk_std/undir_post2wk_mean;

%dir_post4wk_cv = 100*dir_post4wk_std/dir_post4wk_mean;
%undir_post4wk_cv = 100*undir_post4wk_std/undir_post4wk_mean;

dir_post8wk_cv = 100*dir_post8wk_std/dir_post8wk_mean;
undir_post8wk_cv = 100*undir_post8wk_std/undir_post8wk_mean;


%%%%% read normalized FF files and calculate SD %%%%%

norm_pre1 = [pre1(1:length(pre1)-8),'meansubtract.mat'];
norm_pre2 = [pre2(1:length(pre2)-8),'meansubtract.mat'];
norm_post1d = [post1d(1:length(post1d)-8),'meansubtract.mat'];
norm_post2wk = [post2wk(1:length(post2wk)-8),'meansubtract.mat'];
%norm_post4wk = [post4wk(1:length(post4wk)-8),'meansubtract.mat'];
norm_post8wk = [post8wk(1:length(post8wk)-8),'meansubtract.mat'];

% load files for normalized FF
temp = load (norm_pre1);
dir_pre1_norm = getfield(temp,'deltaF_dir'); undir_pre1_norm = getfield(temp,'deltaF_undir'); P_pre1 = getfield(temp,'p'); 
temp = load (norm_pre2);
dir_pre2_norm = getfield(temp,'deltaF_dir'); undir_pre2_norm = getfield(temp,'deltaF_undir'); P_pre2 = getfield(temp,'p'); 
temp = load (norm_post1d);
dir_post1d_norm = getfield(temp,'deltaF_dir'); undir_post1d_norm = getfield(temp,'deltaF_undir'); P_post1d = getfield(temp,'p'); 
temp = load (norm_post2wk);
dir_post2wk_norm = getfield(temp,'deltaF_dir'); undir_post2wk_norm = getfield(temp,'deltaF_undir'); P_post2wk = getfield(temp,'p'); 
%temp = load (norm_post4wk);
%dir_post4wk_norm = getfield(temp,'deltaF_dir'); undir_post4wk_norm = getfield(temp,'deltaF_undir'); P_post4wk = getfield(temp,'p'); 
temp = load (norm_post8wk);
dir_post8wk_norm = getfield(temp,'deltaF_dir'); undir_post8wk_norm = getfield(temp,'deltaF_undir'); P_post8wk = getfield(temp,'p'); 

% calculate SD for normalized FF
dir_pre1_norm_std = std(dir_pre1_norm);
undir_pre1_norm_std = std(undir_pre1_norm);

dir_pre2_norm_std = std(dir_pre2_norm);
undir_pre2_norm_std = std(undir_pre2_norm);

dir_post1d_norm_std = std(dir_post1d_norm);
undir_post1d_norm_std = std(undir_post1d_norm);

dir_post2wk_norm_std = std(dir_post2wk_norm);
undir_post2wk_norm_std = std(undir_post2wk_norm);

%dir_post4wk_norm_std = std(dir_post4wk_norm);
%undir_post4wk_norm_std = std(undir_post4wk_norm);

dir_post8wk_norm_std = std(dir_post8wk_norm);
undir_post8wk_norm_std = std(undir_post8wk_norm);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot mean and SEM of FF
X = [-7 -1 1 14 56]; 
dir_mean = [dir_pre1_mean, dir_pre2_mean, dir_post1d_mean, dir_post2wk_mean, dir_post8wk_mean];
undir_mean = [undir_pre1_mean, undir_pre2_mean, undir_post1d_mean, undir_post2wk_mean, undir_post8wk_mean];
dir_sem = [dir_pre1_sem, dir_pre2_sem, dir_post1d_sem, dir_post2wk_sem, dir_post8wk_sem];
undir_sem = [undir_pre1_sem, undir_pre2_sem, undir_post1d_sem, undir_post2wk_sem, undir_post8wk_sem];

figure; subplot(2,2,1)
errorbar(X, dir_mean, dir_sem, '-or'); hold on
errorbar(X, undir_mean, undir_sem, '-ob')
title(['fundamental frequency of note "',note,'" of ', bird, '  (Red-Dir, Blue-Undir)']) 
axis([-10 60 0.995*min([dir_mean,undir_mean]) 1.005*max([dir_mean,undir_mean])]);
xt =[]; set(gca, 'XTick',xt); 
text(-8, 0.993*min([dir_mean,undir_mean]), 'pre1', 'rotation', -30)
text(-2.5, 0.993*min([dir_mean,undir_mean]), 'pre2', 'rotation', -30)
text(2.5, 0.993*min([dir_mean,undir_mean]), 'post1-2d', 'rotation', -30)
text(14, 0.993*min([dir_mean,undir_mean]), 'post2wk', 'rotation', -30)
%text(4.9,0.993*min([dir_mean,undir_mean]), 'post4wk', 'rotation', -30)
text(56,0.993*min([dir_mean,undir_mean]), 'post8wk', 'rotation', -30)

% plot CV of FF
dir_cv = [dir_pre1_cv, dir_pre2_cv, dir_post1d_cv, dir_post2wk_cv, dir_post8wk_cv];
undir_cv = [undir_pre1_cv, undir_pre2_cv, undir_post1d_cv, undir_post2wk_cv, undir_post8wk_cv];
subplot(2,2,3); plot(X, dir_cv./100, '-or'); hold on
plot(X, undir_cv./100, '-ob')
title(['CV of fundamental frequency'])
Ymax = (1.05*max([dir_cv,undir_cv]))/100
axis([-10 60 0 Ymax]);
xt =[]; set(gca, 'XTick',xt); 
text(-8, Ymax*(-0.05), 'pre1', 'rotation', -30)
text(-2.5, Ymax*(-0.05), 'pre2', 'rotation', -30)
text(2.5, Ymax*(-0.05), 'post1-2d', 'rotation', -30)
text(14, Ymax*(-0.05), 'post2wk', 'rotation', -30)
%text(4.9, Ymax*(-0.05), 'post4wk', 'rotation', -30)
text(56, Ymax*(-0.05), 'post8wk', 'rotation', -30)

% plot SD of FF
dir_std = [dir_pre1_std, dir_pre2_std, dir_post1d_std, dir_post2wk_std, dir_post8wk_std];
undir_std = [undir_pre1_std, undir_pre2_std, undir_post1d_std, undir_post2wk_std, undir_post8wk_std];
subplot(2,2,2); 
plot(X, dir_std, '-or'); hold on
plot(X, undir_std, '-ob')
title(['SD of fundamental frequency'])
Ymax = 1.05*max([dir_std,undir_std])
axis([-10 60 0 Ymax]);
xt =[]; set(gca, 'XTick',xt); 
text(-8, Ymax*(-0.05), 'pre1', 'rotation', -30)
text(-2.5, Ymax*(-0.05), 'pre2', 'rotation', -30)
text(2.5, Ymax*(-0.05), 'post1-2d', 'rotation', -30)
text(14, Ymax*(-0.05), 'post2wk', 'rotation', -30)
%text(4.9, Ymax*(-0.05), 'post4wk', 'rotation', -30)
text(56, Ymax*(-0.05), 'post8wk', 'rotation', -30)

% plot SD of normalized FF
dir_norm_std = [dir_pre1_norm_std, dir_pre2_norm_std, dir_post1d_norm_std, dir_post2wk_norm_std, dir_post8wk_norm_std];
undir_norm_std = [undir_pre1_norm_std, undir_pre2_norm_std, undir_post1d_norm_std, undir_post2wk_norm_std, undir_post8wk_norm_std];
subplot(2,2,4); 
plot(X, dir_norm_std, '-or'); hold on
plot(X, undir_norm_std, '-ob');hold on
title(['SD of % difference from local mean'])
Ymax = 1.05*max([dir_norm_std,undir_norm_std])
axis([-10 60 0 Ymax]);
xt =[]; set(gca, 'XTick',xt); 
text(-8, Ymax*(-0.05), 'pre1', 'rotation', -30)
text(-2.5, Ymax*(-0.05), 'pre2', 'rotation', -30)
text(2.5, Ymax*(-0.05), 'post1-2d', 'rotation', -30)
text(14, Ymax*(-0.05), 'post2wk', 'rotation', -30)
%text(4.9, Ymax*(-0.05), 'post4wk', 'rotation', -30)
text(56, Ymax*(-0.05), 'post8wk', 'rotation', -30)

% plot astarisks for significance
P = [P_pre1, P_pre2, P_post1d, P_post2wk, P_post8wk];
for a = 1:length(P)
    if P(a)<0.05
        plot(X(a), min([dir_norm_std(a),undir_norm_std(a)])-0.1, '*k'); hold on
    end
end

% calculate difference and ratio of normalized SD just in case
norm_std_diff = undir_norm_std - dir_norm_std;
norm_std_ratio = undir_norm_std./dir_norm_std;

% save the figure
fig_name =  ['FF_',note,'_summary_meansubtract.fig']
saveas(gcf,fig_name)

% save the variables
saved_file_name = ['FF_', note, '_summary_meansubtract.mat']
save(saved_file_name,'dir_mean','undir_mean','dir_sem','undir_sem','dir_cv','undir_cv','dir_std','undir_std',...
    'dir_norm_std','undir_norm_std','norm_std_diff','norm_std_ratio','P','bird','note')