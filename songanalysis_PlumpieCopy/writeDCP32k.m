%writess DCP songfile (32kHz) from a vector and name it "AM- ".

function song = writeDCP32k(filename,signal);

SAMP_RATE = 32000;

fid = fopen(filename,'r','b');
header = fread(fid,11,'int16');
fclose(fid);

newfilename = ['AM',filename];
signal = round(signal);
fid = fopen(newfilename,'w','b');
fwrite(fid,header,'int16');
fwrite(fid,signal,'int16');
fclose(fid);
