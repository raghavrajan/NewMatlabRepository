function batch_info = batch_select(batch_info)

if ~isfield(batch_info,'batchfilename') | isempty(batch_info.batchfilename) 
  batch_info.batchfilename = input('Enter batch file name: ','s');
end

while 1
  if isempty(batch_info.batchfilename)
    disp('You must enter a batch file. Try again.')
  else
    batch_info.fid = fopen(batch_info.batchfilename,'r');
    if batch_info.fid == -1
      disp('Invalid batch file name. Try again.') 
    else
      break;
    end
  end
  batch_info.batchfilename = input('Enter batch file name: ','s');
end

if ~isfield(batch_info,'basedir') | isempty(batch_info.basedir) | ~exist(batch_info.basedir,'dir') 
  batch_info.basedir = input('Enter base directory for batch data files [/]: ','s');
end

while 1
  if ~isempty(batch_info.basedir) & ~exist(batch_info.basedir,'dir')
    disp('You must enter a valid directory. Try again.')
  else
    break;
  end
  batch_info.basedir = input('Enter base directory for batch data files [/]: ','s');
end

