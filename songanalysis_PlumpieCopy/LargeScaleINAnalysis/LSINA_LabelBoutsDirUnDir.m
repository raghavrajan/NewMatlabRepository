function [BoutDirUnDir] = LSINA_LabelBoutsDirUnDir(SURecordingDetails)

BoutDirUnDir = [];
% This is a function to go through the recording details and label each
% individual bout as directed (1) or undirected (0)
% If a bout starts before the female went in and ends after the female went
% in, then I will label it as (2) to represent such ambiguous bouts. I will
% use these only if I really need to.

for i = 1:size(SURecordingDetails.Bouts,1),
    % May be easiest to treat all files as continuous irrespective of
    % whether it is triggered or continuous and get one long time base.
    % Then I can just check whether bout onset and offset times fall within
    % the continuous on and off times for directed presentation
    
    BoutOnsetFile = SURecordingDetails.Bouts(i,3);
    BoutOnsetTime = SURecordingDetails.Bouts(i,5);
    
    BoutOffsetFile = SURecordingDetails.Bouts(i,4);
    BoutOffsetTime = SURecordingDetails.Bouts(i,6);
    
    if (BoutOnsetFile > 1)
        ContinuousBoutOnsetTime = sum(SURecordingDetails.FileLen(1:(BoutOnsetFile -1))) + BoutOnsetTime;
        else
        ContinuousBoutOnsetTime = BoutOnsetTime;
    end
        
    if (BoutOffsetFile > 1)
        ContinuousBoutOffsetTime = sum(SURecordingDetails.FileLen(1:(BoutOffsetFile -1))) + BoutOffsetTime;
    else
        ContinuousBoutOffsetTime = BoutOffsetTime;
    end
        
    DirUndir = 0;
    for j = 1:size(SURecordingDetails.DirPresentations,1),
        if (SURecordingDetails.DirPresentations(j,1) > 1)
            ContinuousDirOnsetTime = sum(SURecordingDetails.FileLen(1:(SURecordingDetails.DirPresentations(j,1) - 1))) + SURecordingDetails.DirPresentations(j,5);
        else
            ContinuousDirOnsetTime = SURecordingDetails.DirPresentations(j,5);
        end
        if (SURecordingDetails.DirPresentations(j,4) > 1)
            ContinuousDirOffsetTime = sum(SURecordingDetails.FileLen(1:(SURecordingDetails.DirPresentations(j,4) - 1))) + SURecordingDetails.DirPresentations(j,8);
        else
            ContinuousDirOffsetTime = SURecordingDetails.DirPresentations(j,8);
        end
        
        % Now, to see if the continuous bout onset and offset time fall
        % within this directed presentation or not
        
        if (((ContinuousBoutOnsetTime >= ContinuousDirOnsetTime) & (ContinuousBoutOnsetTime <= ContinuousDirOffsetTime)) || ((ContinuousBoutOffsetTime >= ContinuousDirOnsetTime) & (ContinuousBoutOffsetTime <= ContinuousDirOffsetTime)))
            % Either bout onset of bout offset falls within the directed
            % presentation
            % Now to check whether both fall within, or just one falls
            % within
            if (((ContinuousBoutOnsetTime >= ContinuousDirOnsetTime) & (ContinuousBoutOnsetTime <= ContinuousDirOffsetTime)) && ((ContinuousBoutOffsetTime >= ContinuousDirOnsetTime) & (ContinuousBoutOffsetTime <= ContinuousDirOffsetTime)))
                DirUndir = 1;
            else
                DirUndir = 2; % ambiguous case where only one of bout onset or offset falls within the dir presentation, while the other does not
            end
        end
    end
    BoutDirUnDir(i) = DirUndir;
end
