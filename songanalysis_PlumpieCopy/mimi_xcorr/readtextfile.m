function list = readtextfile(listname);
%reads a text file into a 2-D character array
list = [];
fid = fopen(listname,'r');
ln = 1;
disp('Reading the text file...')

while feof(fid) == 0
   % disp(['Reading line # ', num2str(ln)])
    list = char(list,fgetl(fid));
    ln = ln+1;
end
list(1,:)=[];
fclose(fid);