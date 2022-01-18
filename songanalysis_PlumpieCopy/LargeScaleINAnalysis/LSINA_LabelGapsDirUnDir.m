function [Gaps, GapDirUnDir] = LSINA_LabelGapsDirUnDir(SURecordingDetails)

GapDirUnDir = [];
% This is a function to go through the recording details and label each
% individual gap as directed (1) or undirected (0)
% If a gap starts before the female went in and ends after the female went
% in, then I will label it as (2) to represent such ambiguous gaps. I will
% use these only if I really need to.
Gaps = SURecordingDetails.Gaps;
% For continuous data, sometimes there is a long gap with no singing that
% starts during one dir presentation and then carries on till another dir
% presentation. For such gaps, I will split them into separate gaps based
% on the onset and offset times of dir presentation

for i = 1:size(SURecordingDetails.Gaps,1),
    % May be easiest to treat all files as continuous irrespective of
    % whether it is triggered or continuous and get one long time base.
    % Then I can just check whether gap onset and offset times fall within
    % the continuous on and off times for directed presentation
    
    GapOnsetFile = SURecordingDetails.Gaps(i,4);
    GapOnsetTime = SURecordingDetails.Gaps(i,2);
    
    GapOffsetFile = SURecordingDetails.Gaps(i,5);
    GapOffsetTime = SURecordingDetails.Gaps(i,3);
    
    if (GapOnsetFile > 1)
        ContinuousGapOnsetTime = sum(SURecordingDetails.FileLen(1:(GapOnsetFile -1))) + GapOnsetTime;
        else
        ContinuousGapOnsetTime = GapOnsetTime;
    end
        
    if (GapOffsetFile > 1)
        ContinuousGapOffsetTime = sum(SURecordingDetails.FileLen(1:(GapOffsetFile -1))) + GapOffsetTime;
    else
        ContinuousGapOffsetTime = GapOffsetTime;
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
        
        if (((ContinuousGapOnsetTime >= ContinuousDirOnsetTime) & (ContinuousGapOnsetTime <= ContinuousDirOffsetTime)) || ((ContinuousGapOffsetTime >= ContinuousDirOnsetTime) & (ContinuousGapOffsetTime <= ContinuousDirOffsetTime)))
            % Either bout onset of bout offset falls within the directed
            % presentation
            % Now to check whether both fall within, or just one falls
            % within
            if (((ContinuousGapOnsetTime >= ContinuousDirOnsetTime) & (ContinuousGapOnsetTime <= ContinuousDirOffsetTime)) && ((ContinuousGapOffsetTime >= ContinuousDirOnsetTime) & (ContinuousGapOffsetTime <= ContinuousDirOffsetTime)))
                DirUndir = 1;
            else
                DirUndir = 2; % ambiguous case where only one of bout onset or offset falls within the dir presentation, while the other does not
            end
        end
    end
    GapDirUnDir(i) = DirUndir;
end

% if (SURecordingDetails.Continuousdata == 1)
%     AmbiguousGaps = find(GapDirUnDir == 2);
%     % Either gap onset could be before dir presentation start and gap
%     % offset could be after dir presentation start OR gap onset could be
%     % after dir presentation start and gap offset could be after dir
%     % presentation end
%     for j = 1:size(SURecordingDetails.DirPresentations,1),
%         if (((SURecordingDetails.Gaps(i,6) >= SURecordingDetails.DirPresentations(j,9)) & (SURecordingDetails.Gaps(i,6) <= SURecordingDetails.DirPresentations(j,12))) || ((SURecordingDetails.Gaps(i,7) >= SURecordingDetails.DirPresentations(j,9)) & (SURecordingDetails.Gaps(i,7) <= SURecordingDetails.DirPresentations(j,12))))
%             if (((SURecordingDetails.Gaps(i,6) >= SURecordingDetails.DirPresentations(j,9)) & (SURecordingDetails.Gaps(i,6) <= SURecordingDetails.DirPresentations(j,12))))
%                 % onset is within the dir presentation
%             end
%         end
%     end
% end
