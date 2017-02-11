function xcorr_songenv_batch(batchfile1, batchfile2, birdname);

% read batch files
% modified in 2013; batchfile1 is a file of templates (.mat files)
% batchfile 2 is a file listing song files (.wav or observer)

pre_songs = readtextfile(batchfile1)
post_songs = readtextfile(batchfile2)

% run xcorr_songenv for all combinations
% deblank.m removes whitespace characters from string S;
%   just to clean up the filenames
%
for n=1:length(pre_songs(:,1))  %for each pre_song template
    %presong = [deblank(pre_songs(n,:)),'.DwnSmpl.wav']
    presong=load(pre_songs(n,:));
    
    %template_no=pre_songs(n,:);
    %template_no=template_no(1:length(pre_songs(1,:))-4);
  
    template_file=presong.filename(length(birdname)+2:length(presong.filename));
    display(['pre_song template#',num2str(n), '   from   ', template_file])
    template=presong.motif;
    template_no=n;
    
    %calculate the max xcorr with each post-song bout
    for m=1:length(post_songs(:,1))
        %birdname
       
        postsong = [(post_songs(m,:))]
        %postsong = [deblank(post_songs(m,:)),'.DwnSmpl.wav']
              
        %xcorr_songenv is the function that calculates the cross correlation (mean subtracted)
        [maxcorr_envelop] = xcorr_songenv(template,postsong, birdname,template_file,template_no);
        all_coef((n-1)*length(post_songs(:,1))+m) = maxcorr_envelop
       
    end
end

mean_coef = mean(all_coef)
median_coef=median(all_coef)
std_coef = std(all_coef)
sem_coef = std_coef/(length(all_coef)^(1/2))
%sem_coef = std_coef/sqrt(length(all_coef));

saved_file_name = ['XcorrEnv_summary_',birdname]
save(saved_file_name,'all_coef','mean_coef','std_coef','sem_coef')

