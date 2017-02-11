function [List] = PlotSingingFrequency(BatchFileName)

fid = fopen(BatchFileName,'r');
Songs = [];
Date = [];
List = [];
while (~(feof(fid)))
    tline = fgetl(fid);
    List{(end + 1),1} = tline;
    Cbin_index = strfind(tline,'.cbin');
    RecFilename = [tline(1:Cbin_index),'rec'];
    fid2 = fopen(RecFilename);
    while (~(feof(fid2)))
        tline2 = fgetl(fid2);
        if (strfind(tline2,'File created'))
            Commaindex = strfind(tline2,',');
            Date = [Date; tline2((Commaindex(1) - 3):(Commaindex(end) - 1))];
            Time = tline2((Commaindex(end) + 1):end);
            Colonindex = strfind(Time,':');
            FileTime = str2double(Time(1:(Colonindex(1) - 1))) + (str2double(Time((Colonindex(1)+1):(Colonindex(2) - 1))))/60;
        end
        if (strfind(tline2,'rec end'))
            EqualIndex = strfind(tline2,'=');
            RecordLength = str2double(tline2((EqualIndex + 1):(end - 2)))/1000;
        end
    end
    fclose(fid2);
    Songs(end + 1,:) = [FileTime RecordLength];
    List{end,2} = FileTime;
    if (length(strfind(List{end,1},'_dir')) > 0)
        List{end,3} = 'D';
    else
        List{end,3} = 'U';
    end
end
fclose(fid);

flag = 1;
Total_length = 0;
while (flag)
    FirstDates = strmatch(Date(flag,:),Date,'exact');
    figure;
    set(gcf,'Color','w');
    DirNotes = find([List{FirstDates,3}] == 'D');
    UnDirNotes = find([List{FirstDates,3}] == 'U');
    plot(Songs(FirstDates(DirNotes),1),Songs(FirstDates(DirNotes),2),'r.');
    hold on;
    plot(Songs(FirstDates(UnDirNotes),1),Songs(FirstDates(UnDirNotes),2),'b.');
    set(gca,'FontSize',16);
    title(Date(flag,:),'FontSize',18);
    Total_length = Total_length + length(FirstDates);
    if (Total_length == length(Date))
        flag = 0;
    else
        flag = FirstDates(end) + 1;
    end
end

disp('Finished plotting song times');

