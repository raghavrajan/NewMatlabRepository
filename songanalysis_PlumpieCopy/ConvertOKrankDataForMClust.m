function [] = ConvertOKrankDataForMClust(DirectoryName, FileList, FileType, ChanNo)



number_of_files_to_load = length(filelist);

spike_count = 0;
total_spike_count = 0;
T = [];
wv = [];

for file_number = 1:length(filelist),
    temp = strfind(filelist(file_number).name,'data');
    file_numbers(file_number) = str2num(filelist(file_number).name(temp+4:end));
end

[file_numbers, sorted_file_indices] = sort(file_numbers);
            
filelist = filelist(sorted_file_indices);

for i = 1:number_of_files_to_load,
    filename = filelist(i).name;
    a = nc_getbuffer(filename);
    if (i == 1)
        temp = nc_getall(filename);
        sample_rate = temp.sample_rate.data(first_tetrode_chan + 1)/10000; % have to add 1 to the channel no, as netcdf stores it as 0:63 (the 64 channels), while matlab accesses it as 1:64. My notes are written with reference to the 0:63 mode. the 10000 in the denominator is for MClust which likes timestamps in 100us units
        clear temp;
    end
    spikes_to_get = find(a.id == first_tetrode_chan);
    spikes_to_get_second_chan = find(a.id == (first_tetrode_chan + 1));
    spikes_to_get_third_chan = find(a.id == (first_tetrode_chan + 2));
    spikes_to_get_fourth_chan = find(a.id == (first_tetrode_chan + 3));

    minimum_spikes = min([length(spikes_to_get) length(spikes_to_get_second_chan) length(spikes_to_get_third_chan) length(spikes_to_get_fourth_chan)]);

    if (length(spikes_to_get) > minimum_spikes)
        spikes_to_get(minimum_spikes+1:length(spikes_to_get)) = [];
    end
    count = (total_spike_count + 1):1:(length(spikes_to_get) + total_spike_count);
    
    spike_count = length(spikes_to_get);
    total_spike_count = total_spike_count + length(spikes_to_get);
    
    if (spike_count == 0)
    else
        T(count) = a.timestamp(spikes_to_get);
        wv = [wv ; zeros(length(spikes_to_get),4,64)];
    
        wv(count,1,:) = a.data(spikes_to_get,:);
        wv(count,2,:) = a.data(spikes_to_get+1,:);
        wv(count,3,:) = a.data(spikes_to_get+2,:);
        wv(count,4,:) = a.data(spikes_to_get+3,:);
    end
    
    clear spikes_to_get;
    fprintf('Loaded %i spikes from datafile %s\n',spike_count,filename)
end

fprintf('Loaded %i spikes from %i datafiles\n',total_spike_count,number_of_files_to_load);

outputfile = strcat('analysis_files/',filename,'_chan',int2str(first_tetrode_chan),'.mat');

T = T/sample_rate;
T = T';

save(outputfile,'T','wv');