function [] = ModifyRecFiles(DirectoryName)

% This is a function to change the format line in the rec file
% It changes the first line from 
%format:'okrank .....' to format:'krank .....'

cd(DirectoryName);

if (~exist('OriginalRecFiles','dir'))
    mkdir(DirectoryName, 'OriginalRecFiles');
    movefile('*.rec','OriginalRecFiles');
end


cd('OriginalRecFiles');
RecFiles = dir('*.rec');
cd ../;

for i = 1:length(RecFiles),
    disp(RecFiles(i).name);
    OldRecFilefid = fopen(['OriginalRecFiles/',RecFiles(i).name],'r');
    NewRecFilefid = fopen(RecFiles(i).name,'w');
    while ~(feof(OldRecFilefid))
        tline = fgetl(OldRecFilefid);
        if (strfind(tline,'format: "okrank') > 0)
            fprintf(NewRecFilefid,'format: "krank 20080924"\n');
        else
            fprintf(NewRecFilefid,'%s\n',tline);
        end
    end
    fclose(OldRecFilefid);
    fclose(NewRecFilefid);
end

disp('Done');