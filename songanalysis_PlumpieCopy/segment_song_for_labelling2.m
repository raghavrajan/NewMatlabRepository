function [onsets, offsets]=segment_song_for_labelling2(smooth, Fs, min_int, min_dur, threshold, NoiseMSD)

% segment takes smoothed filtered song and returns vectors of note onsets and offsets
% values are in ms 


%threshold input
notetimes=(smooth==threshold);

%extract index values for note onsets and offsets
h=[1 -1];
temp = zeros(size(notetimes,1), size(notetimes,2));
temp(find(notetimes > 0)) = 1;
trans=conv(h,temp);
onsets=find(trans>0);
offsets=find(trans<0);
if (length(onsets) ~= length(offsets))
    disp('number of note onsets and offsets do not match')
else         
    %eliminate short intervals
    temp_int=(onsets(2:length(onsets))-offsets(1:length(offsets)-1))*1000/Fs;
    real_ints=temp_int>min_int;
    onsets=[onsets(1); nonzeros(onsets(2:length(onsets)).*real_ints)];
    offsets=[nonzeros(offsets(1:length(offsets)-1).*real_ints); offsets(length(offsets))];
    
    %eliminate short notes
    temp_dur=(offsets-onsets)*1000/Fs;
    real_durs=temp_dur>min_dur;
%     small_durs = find(temp_dur<=min_dur);
%        
%     if (length(small_durs) > 0)
%         include_intervals = [];
%         if (small_durs(1) == 1)
%             pre_intervals = (onsets(small_durs(2:end)) - offsets(small_durs(2:end) - 1))*1000/Fs;
%             short_pre_intervals = find(pre_intervals <= 20);
%             for i = 1:length(short_pre_intervals),
%                 include_intervals = [include_intervals (small_durs(short_pre_intervals(i) + 1))];
%             end
%         else
%             pre_intervals = (onsets(small_durs) - offsets(small_durs - 1))*1000/Fs;
%             short_pre_intervals = find(pre_intervals <= 20);
%             for i = 1:length(short_pre_intervals),
%                 include_intervals = [include_intervals (small_durs(short_pre_intervals(i)))];
%             end
%         end
% 
%         if (small_durs(end) == length(onsets))
%             post_intervals = (onsets(small_durs(1:end-1) + 1) - offsets(small_durs(1:end-1)))*1000/Fs;
%             short_post_intervals = find(post_intervals <= 20);
%             for i = 1:length(short_post_intervals),
%                 include_intervals = [include_intervals (small_durs(short_post_intervals(i)))];
%             end
%         else
%             post_intervals = (onsets(small_durs + 1) - offsets(small_durs))*1000/Fs;
%             short_post_intervals = find(post_intervals <= 20);
%             for i = 1:length(short_post_intervals),
%                 include_intervals = [include_intervals (small_durs(short_post_intervals(i)))];
%             end
%         end
%         include_intervals = unique(include_intervals);
% 
%         for i = 1:length(include_intervals)
%             if (include_intervals(i) == 1)
%                 onsets(include_intervals(i) + 1) = onsets(include_intervals(i));
%             else
%                 if (include_intervals(i) == length(onsets))
%                     offsets(include_intervals(i) - 1) = offsets(include_intervals(i));
%                 else
%                     pre_interval = onsets(include_intervals(i)) - offsets(include_intervals(i) - 1);
%                     post_interval = onsets(include_intervals(i) + 1) - offsets(include_intervals(i));
%                     if (pre_interval <= post_interval)
%                         offsets(include_intervals(i) - 1) = offsets(include_intervals(i));
%                     else
%                         onsets(include_intervals(i) + 1) = onsets(include_intervals(i));
%                     end
%                 end
%             end
%         end
%     end
    
    onsets=[nonzeros((onsets).*real_durs)];
    offsets=[nonzeros((offsets).*real_durs)];
    
    %convert to ms: peculiarities here are to prevent rounding problem
    % if t_ons is simply replaced with onsets, everything gets rounded
    onsets = onsets*1000/Fs;
    offsets = offsets*1000/Fs;
end

% T = (1:1:length(smooth)) * 1000/Fs;
% MinusSD = T(find(smooth <= (NoiseMSD(1) + (ThreshMultiplier*NoiseMSD(2)))));
% for i = 1:length(onsets),
%     RealOnset = MinusSD(find(MinusSD <= onsets(i), 1, 'last'));
%     if (length(RealOnset > 0))
%         onsets(i) = RealOnset;
%     end
%     RealOffset = MinusSD(find(MinusSD >= offsets(i), 1, 'first'));
%     if (length(RealOffset > 0))
%         offsets(i) = RealOffset;
%     end
% end
[onsets, OnsetIndices] = unique(onsets);
offsets = offsets(OnsetIndices);

[offsets, OffsetIndices] = unique(offsets);
onsets = onsets(OffsetIndices);
    