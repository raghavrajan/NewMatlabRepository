function [] = plot_song_raster_pst(pathname,filename,file_numbers)

figure;
hold on;
raster_y_value = 1;
for i = 1:length(file_numbers),
    song_note_file = strcat(filename,'.',num2str(file_numbers(i)),'.cbin.not.mat');
    spike_time_file = strcat(filename,'.',num2str(file_numbers(i)),'.cbin.spiketimes');
    
    song_notes = load(song_note_file);
    spike_times = load(spike_time_file);
    songlabels = song_notes.labels
    
    bout_beginning = strfind(songlabels,'a');
    bout_end = strfind(songlabels,'e');
    
    if (length(bout_beginning) == length(bout_end))
        for bout_no = 1:length(bout_end),
            bout = [song_notes.onsets(bout_beginning(bout_no)) song_notes.offsets(bout_end(bout_no))];
            trial_spike_times = spike_times(find((spike_times > (bout(1)- 200)) & (spike_times < (bout(2) + 200))));
            trial_spike_times = trial_spike_times - bout(1);
            raster_y_values = ones(length(trial_spike_times),1)*raster_y_value;
            plot(trial_spike_times,raster_y_values,'k+');
            raster_y_value = raster_y_value + 1;
        end
    else
        if (length(bout_beginning) > length(bout_end))
            test = length(bout_beginning) - length(bout_end);
            while (test ~= 0)
                for bout_no = 1:length(bout_end),
                    temp = find((bout_beginning > bout_beginning(bout_no)) & (bout_beginning < bout_end(bout_no)));
                    if (length(temp) > 0)
                        break;
                    end
                end
                bout_beginning(temp) = [];
                test = length(bout_beginning) - length(bout_end);
            end
            for bout_no = 1:length(bout_end),
                bout = [song_notes.onsets(bout_beginning(bout_no)) song_notes.offsets(bout_end(bout_no))];
                trial_spike_times = spike_times(find((spike_times > (bout(1)- 200)) & (spike_times < (bout(2) + 200))));
                trial_spike_times = trial_spike_times - bout(1);
                raster_y_values = ones(length(trial_spike_times),1)*raster_y_value;
                plot(trial_spike_times,raster_y_values,'k+');
                raster_y_value = raster_y_value + 1;
            end
        else
            test = length(bout_beginning) - length(bout_end);
            while (test ~= 0)
                for bout_no = 1:length(bout_beginning),
                    temp = find((bout_end > bout_beginning(bout_no)));
                    if (length(temp) == 0)
                        break;
                    end
                end
                bout_beginning(temp) = [];
                test = length(bout_beginning) - length(bout_end);
            end
            for bout_no = 1:length(bout_end),
                bout = [song_notes.onsets(bout_beginning(bout_no)) song_notes.offsets(bout_end(bout_no))];
                trial_spike_times = spike_times(find((spike_times > (bout(1)- 200)) & (spike_times < (bout(2) + 200))));
                trial_spike_times = trial_spike_times - bout(1);
                raster_y_values = ones(length(trial_spike_times),1)*raster_y_value;
                plot(trial_spike_times,raster_y_values,'k+');
                raster_y_value = raster_y_value + 1;
            end
        end
    end
end    
    