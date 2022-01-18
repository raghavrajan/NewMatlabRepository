function [NewStatusString] = ReactionTimeUpdateStatusString(StatusString, StatusText)

NewStatusString{1} = StatusText;
NewStatusString = [NewStatusString; StatusString];

if (length(NewStatusString) > 35)
    NewStatusString(36:end) = [];
end