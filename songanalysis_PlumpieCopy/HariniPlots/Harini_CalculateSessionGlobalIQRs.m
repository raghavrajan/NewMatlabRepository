function [SessionIQR, Global IQR] = Harini_CalculateSessionGlobalIQRs(Session, MeasureType)

GlobalIQR = [];
SessionIQR = [];

for i = 1:length(Session),
    for j = 1:length(Syll),
        GlobalIQR(end+1) = iqr(eval(['Session(i).Syll(j).Syll', MeasureType]));
        TempSessionIQR = [];
        for RecordingDay = unique(Session(i).Syll(j).RecordingDayIndex),
            for Condition = unique(Session(i).Syll(j).ConditionIndex),
                SessionSylls = find((Session(i).Syll(j).RecordingDayIndex == RecordingDay) & (Session(i).Syll(j).ConditionIndex == ConditionIndex));
                if (~isempty(SessionSylls))
                    