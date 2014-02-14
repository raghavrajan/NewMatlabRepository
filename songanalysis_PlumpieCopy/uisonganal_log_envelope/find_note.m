function find_note(note,expr)


%until the end of batchfile, find and return songs which contain note that matches expression
%note is a string label for a note
%expr is a string containing an expression that can be interpreted by matlab find command (e.g. '>60') 


%get name of metafile containing notefile names

meta_fid = -1;
metafile = 0;
while meta_fid == -1 | metafile == 0 | isempty(metafile)
   disp('select batchfile');
   [metafile, pathname]=uigetfile('*','select batchfile')
   meta_fid=fopen([pathname, metafile]);
   if meta_fid == -1 | metafile == 0
      disp('cannot open file' )
      disp (metafile)
   end
end


while 1
   %get notefile name
     notefile = fscanf(meta_fid,'%s',1);
   %end when there are no more notefiles 
     if (isempty(notefile))
        disp('End of notefiles')
        break
     end
   
   %if notefile exists, get it
     if exist([pathname, notefile])   
       load([pathname, notefile]);
     else
       disp(['cannot find ', notefile])
     end
     
     
     %make sure labels is a column vector
   [x,y]=size(labels);
   if y > x
      labels=labels';
   end 
     
   %check to see if this file contains a note of specified type that satisfies expr
     durations =  get_durs(onsets,offsets,labels,note);
     test = eval(['find(durations',expr,')']);
     if nnz(test > 0)
        %disp(notefile)
        disp(['In ',notefile,': ',note,' ',expr,' at occurance # ',mat2str(test)])
        disp(['   ',labels'])
     end
 
end
     
