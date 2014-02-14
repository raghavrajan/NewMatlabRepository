function ffscan(batch_in)

%quick way to scan through a batch of songs to look at the ff across trials
%requires a batch file w/ a list of the song files and the sampling rate
%Fs=32000 by default

fs=input('Sampling rate?   (default=32kHz)   ');
if length(fs)==0
    fs=32000
end;
fs

%open batch_file
meta_fid=fopen([batch_in]);
if meta_fid == -1 | batch_in == 0
      disp('cannot open file' )
      disp (batch_in)
      return
end

h=figure;

while 1
   
    songfile=fscanf(meta_fid,'%s',1);
    %end when there are no more song files
       if isempty(songfile);
           break
       end
       
     %get the songfile if it exists
     if exist([songfile]);  
        disp(songfile);
       % [filtsong,Fs]=read_filt(songfile);
       filtsong = wavread(songfile);
     else
       disp(['cannot find ', songfile])
     end
     
     
     %plot the specgram and the ff estimates
    
     ffplot(filtsong,fs)
     
     pause  %wait for user to hit any key
     disp(['Press any key to continue'])
     
 end   %while

