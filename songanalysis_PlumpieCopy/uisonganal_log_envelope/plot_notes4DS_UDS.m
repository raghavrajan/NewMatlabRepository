function plot_notes4DS_UDS(notes, elts, dur_data)

%plots up raw durations for the notes listed in the string "notes"
%elts contains the labels for each row of durations in dur_data

%elts and dur_data must have same number of rows
if length(elts) ~= size(dur_data, 1)
   disp('number of note labels must = number of data rows')
   return
end   


%deal with version compatibility issues...
ver_num=version('-release');
ver_num=str2num(ver_num);
if ver_num >= 13
    xticks='xticklabel';
else
    xticks='xticklabels';
end


%open and hold figure
h_durs = figure;
hold on
set(gca,xticks,[])

for i = 1:length(notes)
   %find row number of data labeled by current note
     note_row = find(elts == notes(i));
   %find durations of current note, number of this note, mean and sdev
     durs = dur_data(note_row,:);
     durs = nonzeros(durs);
     n = length(durs);
     %if there were no notes, skip data display
     if n ~= 0
        mean_dur = mean(durs);
        std_dur = std(durs);
       %display data to screen
        disp([notes(i),'   ',num2str(mean_dur),'   ',num2str(std_dur),'   ',num2str(n)]);
       %plot data
        plot(i*ones(size(durs)),durs,'+r')
        line([i-.1,i+.1],[mean_dur, mean_dur],'color','g')
        line([i-.1,i+.1],[mean_dur+std_dur, mean_dur+std_dur],...
             'color','g',...
             'linestyle',':')
        line([i-.1,i+.1],[mean_dur-std_dur, mean_dur-std_dur],...
             'color','g',...
             'linestyle',':')
        text(i+0.02,min(durs),['(', num2str(n), ')'],...
             'horizontalalignment','center')
        ylim = get(gca,'ylim');
        ylim(2) = max(ylim(2), max(durs)+8); 
        set(gca,'ylim',[0 ylim(2)]);
     end
        
     %set label on x axis
     set(gca,'xtick',[1:i])
     old_xticklabels = get(gca,xticks);
     set(gca,xticks,[old_xticklabels; notes(i)])
                            
end

hold off

