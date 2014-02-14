function xcorr_songenv_batch(batchfile1, batchfile2, birdname);

% read batch files
pre_songs = readtextfile(batchfile1)
post_songs = readtextfile(batchfile2)

% run xcorr_songenv for all combinations
for n=1:length(pre_songs(:,1))
    presong = [deblank(pre_songs(n,:)),'.DwnSmpl.wav']
    for m=1:length(post_songs(:,1))
        postsong = [deblank(post_songs(m,:)),'.DwnSmpl.wav']
        birdname
        [maxcorr_envelop] = xcorr_songenv(presong,postsong, birdname);
        all_coef((n-1)*length(post_songs(:,1))+m) = maxcorr_envelop;
    end
end

mean_coef = mean(all_coef)
std_coef = std(all_coef)
sem_coef = std_coef/(length(all_coef)^(1/2))

saved_file_name = ['XcorrEnv_summary_',birdname]
save(saved_file_name,'all_coef','mean_coef','std_coef','sem_coef')

