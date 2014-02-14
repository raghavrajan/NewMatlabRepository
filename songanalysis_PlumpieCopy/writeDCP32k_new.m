%writess DCP songfile (32kHz) from a vector and name it a prefix.

function song = writeDCP32k_new(filename,signal,prefix);

SAMP_RATE = 32000;

fid = fopen(filename,'r','b');
header = fread(fid,11,'int16');
fclose(fid);

newfilename = [prefix,filename];
signal = round(signal);
fid = fopen(newfilename,'w','b');
fwrite(fid,header,'int16');
fwrite(fid,signal,'int16');
fclose(fid);
