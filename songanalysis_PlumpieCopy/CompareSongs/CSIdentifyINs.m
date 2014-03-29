function [INs] = CSIdentifyINs(Labels, MotifSylls, INLabels, MotifInitiationSylls)

MotifStarts = [];
for i = 1:length(MotifInitiationSylls),
    MotifStarts = [MotifStarts; find(Labels == MotifInitiationSylls(i))];
end

MotifStarts = sort(MotifStarts);
MotifStarts(find(diff(MotifStarts) < 2) + 1) = [];

if (~isempty(MotifStarts))
    INs.BeginningSyllNo = MotifStarts(1) - 1;
    INs.BeginningSylls = 1:1:MotifStarts(1) - 1;
    INs.BeginningSyllPos = -(MotifStarts(1) - 1):1:-1;
    
    INs.INsBeginningINs = [];
    INs.INsPosition = [];
    for i = (MotifStarts(1) - 1):-1:1,
        if (~isempty(find(INLabels == Labels(i))))
           INs.INsBeginningINs = [INs.INsBeginningINs; i];
        else
            break;
        end
    end
    INs.INsBeginningINs = sort(INs.INsBeginningINs);
    INs.INsBeginningINPos = -(INs.INsBeginningINs(end) - INs.INsBeginningINs(1) + 1):1:-1;
    INs.INsBeginningINPos = -(INs.INsBeginningINs(end) - INs.INsBeginningINs(1) + 1):1:-1;
    INs.INsBeginningINNo = length(INs.INsBeginningINs);
else
    INs = [];
end
        
