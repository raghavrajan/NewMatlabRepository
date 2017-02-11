function cbs_time_extract(Dir_Batch, Undir_Batch); 

% read files for recording time
dir_time = readtextfile('time_dir');
n_dir_time = length(dir_time(:,1));

undir_time = readtextfile('time_undir');
n_undir_time  = length(undir_time(:,1));

% read batch files
dir_batch = readtextfile(Dir_Batch);
n_dir_batch = length(dir_batch(:,1));

undir_batch = readtextfile(Undir_Batch);
n_undir_batch = length(undir_batch(:,1));


% extract a recording time of the dir_batch file from the dir_time file.
dir_time_batch = cell(n_dir_batch,1);
i = 1;
for n=1:n_dir_batch
    filename = deblank(dir_batch(n,:));
    for m=1:n_dir_time
        loc_filename = findstr(dir_time(m,:), filename);
        if loc_filename>0
           line1= deblank(dir_time(m,:));
           loc_Duration = findstr(line1, 'Duration');
           dir_time_batch(i) = {[line1(loc_Duration-15:loc_Duration-11),' ',filename]}
        end
    end
    i = i+1;
end

% extract a recording time of the undir_batch file from the undir_time file.
undir_time_batch = cell(n_undir_batch,1);
i = 1;
for n=1:n_undir_batch
    filename = deblank(undir_batch(n,:));
    for m=1:n_undir_time
        loc_filename = findstr(undir_time(m,:), filename);
        if loc_filename>0
           line1= deblank(undir_time(m,:));
           loc_Duration = findstr(line1, 'Duration');
           undir_time_batch(i) = {[line1(loc_Duration-15:loc_Duration-11),' ',filename]}
        end
    end
    i = i+1;
end

% sort the cbs file names by recording time and save them into time_combine.batch.
comb_time_batch = [dir_time_batch; undir_time_batch];
comb_time_batch2 = sort(comb_time_batch)

fid = fopen ('time_combine.batch','w');
for p=1:length(comb_time_batch2)
    fprintf(fid,'%s\n',comb_time_batch2{p});
end
fclose(fid);
