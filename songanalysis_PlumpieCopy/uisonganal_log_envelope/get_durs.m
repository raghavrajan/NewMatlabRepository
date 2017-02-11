function [durations] = get_durs(onsets, offsets, labels, notes)

%returns durations for the notes in the input string 'notes'
%if no note string is provided, returns all durations
 

durations = [];

%if no notes are provided returns all durations
if nargin <= 3
  durations = offsets - onsets; 
else
   for i = 1:length(notes)
    if ~isempty(labels)  
     index = find(notes(i) == labels);
     new_durs = offsets(index) - onsets(index); 
     durations = [durations; new_durs];
    end 
  end
end  
