function batch_info = batch_read(batch_info)

ifile = 0;
while 1
  filename = fgetl(batch_info.fid);
  % Check for whitespace here!
  spaceflag = 0;
  if isspace(filename)
    spaceflag = 1
  end
  if (~ischar(filename)) | isempty(filename) | spaceflag
    disp('End of batch file reached.')
    break
  end
  ifile = ifile+1;  
  batch_info.filenames{ifile} = filename;
end

batch_info.nfiles = ifile;

fclose(batch_info.fid);
