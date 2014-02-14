function [NewStatusString] = UpdateStatusString(StatusString, StatusText)

if (size(StatusString, 1) == 1)
    NewStatusString{1} = StatusString;
else
    NewStatusString = StatusString;
end

if (iscell(StatusText))
    for i = 1:length(StatusText),
        NewStatusString{end + 1} = StatusText{i};
    end
else
    NewStatusString{end + 1} = StatusText;
end

if (length(NewStatusString) > 35)
    Difference = (length(NewStatusString) - 35);
    NewStatusString(2:(2 + Difference - 1)) = [];
end