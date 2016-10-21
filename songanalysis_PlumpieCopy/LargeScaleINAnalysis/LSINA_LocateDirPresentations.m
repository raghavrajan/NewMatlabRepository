function [DirPresentations] = LSINA_LocateDirPresentations(SURecordingDetails)

DirPresentations = [];

% Now return a variable that has one row for each directed presenation.
% The important events are door open(female in), door close, door
% open(female out) and door close. The first 4 columns will have the file
% nos. in which these events occur and the next 4 columns will have the
% times within that file when the events occur.

for i = 1:length(SURecordingDetails.NoteInfo),
    % See if there are any female in or female out or door close events in
    % the labels for each file
    FemaleInEvents = find(SURecordingDetails.NoteInfo{i}.labels == SURecordingDetails.FemaleinLabel);
    if (isempty(FemaleInEvents))
        FemaleIn(i) = 0;
    else
        FemaleIn(i) = FemaleInEvents;
    end
    
    FemaleOutEvents = find(SURecordingDetails.NoteInfo{i}.labels == SURecordingDetails.FemaleoutLabel);
    if (isempty(FemaleOutEvents))
        FemaleOut(i) = 0;
    else
        FemaleOut(i) = FemaleOutEvents;
    end
    
    CloseDoorEvents = find(SURecordingDetails.NoteInfo{i}.labels == SURecordingDetails.Closedoorlabel);
    if (isempty(CloseDoorEvents))
        CloseDoor(i) = 0;
    else
        CloseDoor(i) = CloseDoorEvents;
    end
end    

% Now I will look for female in first assuming that is there and will then
% check that there are no female outs before that
FemaleInEvents = find(FemaleIn);
FemaleOutEvents = find(FemaleOut);
CloseDoorEvents = find(CloseDoor);

for i = 1:length(FemaleInEvents),
    % First check if there are any female out events before the first
    % female in event
    if (i == 1)
        FemaleOutBeforeFirstFemaleIn = find(FemaleOutEvents < FemaleInEvents(i));
        CloseDoorBeforeFirstFemaleIn = find(CloseDoorEvents < FemaleInEvents(i));
        
        if (~isempty(FemaleOutBeforeFirstFemaleIn))
            
            % If there are then assume that the female came in at time 0 in the
            % first file - door close also at the same time
            % First put indices of files
            DirPresentations(end+1,1:4) = [1 1 FemaleOutEvents(FemaleOutBeforeFirstFemaleIn) CloseDoorEvents(CloseDoorBeforeFirstFemaleIn)];
            % Now put times of events
            FemaleOutTime = SURecordingDetails.NoteInfo{FemaleOutEvents(FemaleOutBeforeFirstFemaleIn)}.onsets(FemaleOut(FemaleOutEvents(FemaleOutBeforeFirstFemaleIn)));
            CloseDoorTime = SURecordingDetails.NoteInfo{CloseDoorEvents(CloseDoorBeforeFirstFemaleIn)}.onsets(CloseDoor(CloseDoorEvents(CloseDoorBeforeFirstFemaleIn)));
            DirPresentations(end,5:8) = [0 0 FemaleOutTime CloseDoorTime];
        end
    end
    
    NextFemaleOut = find(FemaleOutEvents >= FemaleInEvents(i), 1, 'first');
    NextCloseDoorEvent = find(CloseDoorEvents >= FemaleInEvents(i), 1, 'first');
    
    FemaleInTime = SURecordingDetails.NoteInfo{FemaleInEvents(i)}.onsets(FemaleIn(FemaleInEvents(i)));
    CloseDoorTime = SURecordingDetails.NoteInfo{CloseDoorEvents(NextCloseDoorEvent)}.onsets(CloseDoor(CloseDoorEvents(NextCloseDoorEvent)));
        
    if (isempty(NextFemaleOut))
        % Then there is no female out and I will assume that female was
        % taken out at the end of the last file
        DirPresentations(end+1,1:4) = [FemaleInEvents(i) CloseDoorEvents(NextCloseDoorEvent) length(SURecordingDetails.NoteInfo) length(SURecordingDetails.NoteInfo)];
        % Now to put times of events
        DirPresentations(end,5:8) = [FemaleInTime CloseDoorTime SURecordingDetails.FileLen(end) SURecordingDetails.FileLen(end)];
    else
        FemaleOutCloseDoorEvent = find(CloseDoorEvents >= FemaleOutEvents(NextFemaleOut), 1, 'first');
        DirPresentations(end+1,1:4) = [FemaleInEvents(i) CloseDoorEvents(NextCloseDoorEvent) FemaleOutEvents(NextFemaleOut) CloseDoorEvents(FemaleOutCloseDoorEvent)];
        % Now to put times of events
        FemaleOutTime = SURecordingDetails.NoteInfo{FemaleOutEvents(NextFemaleOut)}.onsets(FemaleOut(FemaleOutEvents(NextFemaleOut)));
        FemaleOutCloseDoorTime = SURecordingDetails.NoteInfo{CloseDoorEvents(FemaleOutCloseDoorEvent)}.onsets(CloseDoor(CloseDoorEvents(FemaleOutCloseDoorEvent)));
        DirPresentations(end,5:8) = [FemaleInTime CloseDoorTime FemaleOutTime FemaleOutCloseDoorTime];
    end
end
disp('Finished');