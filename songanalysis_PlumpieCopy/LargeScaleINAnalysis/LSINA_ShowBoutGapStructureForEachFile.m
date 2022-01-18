function [] = LSINA_ShowBoutGapStructureForEachFile(SURecordingDetails)

Index = 1;
figure;
hold on;
BoutColours = 'bry';
GapColours = 'cmg';

GapYPos = 0.2;
BoutYPos = 0.2;

set(gca, 'DefaultLineLineWidth', 4);


% First plot bouts at -0.1 below the file #
for i = 1:size(SURecordingDetails.Bouts,1)
    if (SURecordingDetails.Bouts(i,3) == SURecordingDetails.Bouts(i,4))
        plot([SURecordingDetails.Bouts(i,5) SURecordingDetails.Bouts(i,6)], [SURecordingDetails.Bouts(i,3) SURecordingDetails.Bouts(i,4)] - BoutYPos, BoutColours(SURecordingDetails.BoutDirUnDir(i) + 1));
    else
        for j = SURecordingDetails.Bouts(i,3):SURecordingDetails.Bouts(i,4),
            if (j == SURecordingDetails.Bouts(i,3))
                plot([SURecordingDetails.Bouts(i,5) SURecordingDetails.FileLen(j)], [j j] - BoutYPos, BoutColours(SURecordingDetails.BoutDirUnDir(i) + 1));
            else
                if (j == SURecordingDetails.Bouts(i,4))
                    plot([0 SURecordingDetails.Bouts(i,6)], [j j] - BoutYPos, BoutColours(SURecordingDetails.BoutDirUnDir(i) + 1));
                else
                    plot([0 SURecordingDetails.FileLen(j)], [j j] - BoutYPos, BoutColours(SURecordingDetails.BoutDirUnDir(i) + 1));
                end
            end
        end
    end
end

% Next plot gaps at +0.1 above file #
for i = 1:size(SURecordingDetails.Gaps,1)
    if (SURecordingDetails.Gaps(i,4) == SURecordingDetails.Gaps(i,5))
        plot([SURecordingDetails.Gaps(i,2) SURecordingDetails.Gaps(i,3)], [SURecordingDetails.Gaps(i,4) SURecordingDetails.Gaps(i,5)] + GapYPos, GapColours(SURecordingDetails.GapDirUnDir(i) + 1));
    else
        for j = SURecordingDetails.Gaps(i,4):SURecordingDetails.Gaps(i,5),
            switch j
                case SURecordingDetails.Gaps(i,4)
                    plot([SURecordingDetails.Gaps(i,2) SURecordingDetails.FileLen(j)], [j j] + GapYPos, GapColours(SURecordingDetails.GapDirUnDir(i) + 1));

                case SURecordingDetails.Gaps(i,5)
                    plot([0 SURecordingDetails.Gaps(i,3)], [j j] + GapYPos, GapColours(SURecordingDetails.GapDirUnDir(i) + 1));
                    
                otherwise
                    plot([0 SURecordingDetails.FileLen(j)], [j j] + GapYPos, GapColours(SURecordingDetails.GapDirUnDir(i) + 1));
            end
        end
    end
end

% Finally plot the directed presentation times at the file #
for i = 1:size(SURecordingDetails.DirPresentations,1),
    if (SURecordingDetails.DirPresentations(i,1) == SURecordingDetails.DirPresentations(i,4))
        plot([SURecordingDetails.DirPresentations(i,5) SURecordingDetails.DirPresentations(i,8)], [SURecordingDetails.DirPresentations(i,1) SURecordingDetails.DirPresentations(i,4)], 'r');
    else
        for j = SURecordingDetails.DirPresentations(i,1):SURecordingDetails.DirPresentations(i,4),
            switch j
                case SURecordingDetails.DirPresentations(i,1)
                    plot([SURecordingDetails.DirPresentations(i,5) SURecordingDetails.FileLen(j)], [j j], 'k');
                    
                case SURecordingDetails.DirPresentations(i,4)
                    plot([0 SURecordingDetails.DirPresentations(i,8)], [j j], 'k');
                    
                otherwise
                    plot([0 SURecordingDetails.FileLen(j)], [j j], 'k');
            end
        end
    end
end

axis tight;
xlabel('Time (sec)');
ylabel('File #');
Temp = axis;
set(gca, 'XTick', [0:3000:Temp(2)], 'XTickLabel', 0:3:Temp(2)/1000);