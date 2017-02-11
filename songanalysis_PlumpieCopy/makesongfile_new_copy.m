function makesongfile_new_copy(infile,outfile);

%makesongfile(infile,outfile);
%
%	MAKESONGFILE dumps vector INFILE to a observer-
%	playable file OUTFILE in directory d:/songs/!!!
%  (OUTFILE should be name in single quotes)ma
%

%infile=(infile/max(infile))*100;  %may not be necc!
filename=['/home/raghav/',outfile];


[fp,m1]=fopen(filename,'w','b');
fwrite(fp,infile,'int16');
fclose(fp);
