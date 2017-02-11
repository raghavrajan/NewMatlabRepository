% This was given by Brian, and modified by Patrick.
% This converts a wav file into a dcp file and attenuates the amplitude (-10dB)

MAXSHORT = 32767
MINSHORT = -32768
MIN18BIT = -131072 
MAX18BIT = 131071

%songfile = input('Enter name of song wav file: ', 's');
songtype = input(['Enter output song file type: (1) dcp 32kHz 16 bit, ' ...
      ' (2) dcp 40kHz 32 bit: '])
%songwfile = input('Enter name of output song wav file: ', 's');

if songtype == 1
  nwbits = 16
  wtype = 'int16'
  Fsw = 32000
elseif songtype == 2
  nwbits = 32;
  wtype = 'int32'
  Fsw = 40000
end

% Will assume cbs files recorded at this rate!
Fs = 44100

% Maximum length of dcp song in sec
MAXSONGDUR = 8.0

% Resampling algorithm used.
resamptype = 'matlab'

nfft = 2048;

% Normalization parameters
cutfrac = 0.025
Nstd = 5.0

disp('Normalization Types Menu:')
disp('(0) [None]')
disp('(1) Center cut peak RMS')
disp('(2) Center cut RMS')
disp('(3) RMS')
disp('(4) Set RMS')
disp('(5) Set Range')
disp('(6) MAX Range')
disp(' ')
normcode = input(['What type of normalization do you want: (0-6)? '])
if normcode == 1
  NORMTYPE = 'centercut_peakrms'
elseif normcode == 2
  NORMTYPE = 'centercut_rms'
elseif normcode == 3
  NORMTYPE = 'rms'
elseif normcode == 4
  NORMTYPE = 'set_rms'
  setrmsval = input('Enter rms value for normalized data: ')
elseif normcode == 5
  NORMTYPE = 'set_range'
  setrangeval = input('Enter positive range for normalized data: ')
elseif normcode == 6
  NORMTYPE = 'max'
else
  NORMTYPE = 'none'
end

% Filtering parameters
disp('Bandpass Filter Types Menu:')
disp('(0) [None]')
disp('(1) Kaiser window FIR (500-8000 Hz)')
disp('(2) Hanning Window FIR (variable Fc)')
disp('(3) Butterworth IIR (variable Fc)')
disp('(4) High pass Hanning window FIR (variable Fc)')
disp(' ')
filtcode = input(['What type of filter do you want: (0-4)? '])
if filtcode == 1
  FILTTYPE = 'kaiserfir'
elseif filtcode == 2
  FILTTYPE = 'hanningfir'
  F_low = input('Enter low frequency cutoff: ')
  F_high = input('Enter high frequency cutoff: ')  
elseif filtcode == 3
  FILTTYPE = 'butter'
  F_low = input('Enter low frequency cutoff: ')
  F_high = input('Enter high frequency cutoff: ')  
elseif filtcode == 4
  FILTTYPE = 'hipass'
  F_low = input('Enter low frequency cutoff: ')
else
  FILTTYPE = 'none'
end
  
if filtcode
  yesno = input('Do you want to use matlab ''filtfilt'' function (''filter'' is default) (y/N)? ','s')
  if strncmpi(yesno,'y',1)
    do_filtfilt = 1
  else
    do_filtfilt = 0
  end
  switch lower(FILTTYPE)
    case 'kaiserfir'
      % FIR filter design with a Kaiser window
      % This one looks like a nice filter to use (don't erase!)
      Rp = 3.266
      Rs = 30
      %fbands = [500 600 8000 8800]
      fbands = [293 453.1 8223 8328]
      amps = [0 1 0]
      devs = [10^(-Rs/20) (10^(Rp/20)-1)/(10^(Rp/20)+1) 10^(-Rs/20)]
      [nfir,Wnfir,beta,ftype] = kaiserord(fbands,amps,devs,Fs);
      nfir = nfir + rem(nfir,2)
      ndelay = fix(nfir/2)
      bfir = fir1(nfir,Wnfir,ftype,kaiser(nfir+1,beta));
      %figure;
      %freqz(bfir,1,nfft,Fs)
      %figure;
      %grpdelay(bfir,1,nfft,Fs)
    case 'hanningfir'
      nfir = 512
      ndelay = fix(nfir/2)
      bfir = fir1(nfir,[F_low*2/Fs, F_high*2/Fs]);
      %figure;
      %freqz(bfir,1,nfft,Fs)
      %figure;
      %grpdelay(bfir,1,nfft,Fs)
   case 'hipass'
      nfir = 512
      ndelay = fix(nfir/2) 
      bfir = fir1(nfir, F_low*2/Fs, 'high');
      %figure;
      %freqz(bfir,1,nfft,Fs)
      %figure;
      %grpdelay(bfir,1,nfft,Fs)      
   end
end

% Get batch file information

while 1
  batchfilename = input('Enter batch file name: ','s');
  if isempty(batchfilename)
    disp('You must enter a batch file. Try again.')
  else
    fid = fopen(batchfilename,'r');
    if fid == -1
      disp('Invalid batch file name. Try again.') 
    else
      break;
    end
  end
end
songbasedir = input('Enter base directory for song files: ','s')

ifile = 0
while 1
  songname = fgetl(fid);
  % Check for whitespace here!
  spaceflag = 0;
  if isspace(songname)
    spaceflag = 1
  end
  if (~ischar(songname)) | isempty(songname) | spaceflag
    disp('End of batch file reached.')
    break
  end
  songfile = fullfile(songbasedir,songname);
  [song,Fs,nbits]=wavread(songfile);
  songlen = length(song)
  
  %  fids = fopen(songfile,'r','b');
  %  [song, songlen] = fread(fids,inf,'int16');
  %  fclose(fids);
  ifile = ifile+1  

  if filtcode
    switch lower(FILTTYPE)
      case 'butter'
	filtsong = bandpass(song,Fs,F_low,F_high);  
      case {'hanningfir', 'kaiserfir', 'hipass'}
	if do_filtfilt
	  filtsong = filtfilt(bfir,1,song);
	else
	  z = [];
	  [filtsongd, z] = filter(bfir,1,song,z);
	  % Note filtered song corrected for the delay is shorter in length!
	  filtsong = filtsongd(1+ndelay:length(filtsongd));
	end
    end
  else 
    filtsong = song;    
  end

  % Normalize song to maximize dynamic range
  % Three way to try to do this
  
  songmax = max(abs(filtsong))
  songmean = mean(filtsong)
  songstd = std(filtsong)
  disp(['Song mean: ', num2str(songmean), ' std: ', num2str(songstd)])
  normsong = filtsong - songmean;
  ccut = songmax*cutfrac
  normsonglen = length(filtsong);

  switch lower(NORMTYPE)
    case 'centercut_peakrms'
      % This is what FET's songfilt used to do
      ccrms = 0;
      ncclip = 0;
      for i=2:normsonglen-1
	if (normsong(i) > ccut & normsong(i) > normsong(i-1) & ...
	      normsong(i) > normsong(i+1)) | ...
	      (normsong(i) < -ccut & normsong(i) < normsong(i-1) & ...
	      normsong(i) < normsong(i+1)) 
	  ccrms = ccrms + normsong(i)^2;
	  ncclip = ncclip + 1;
	end  
      end
      ccrms = sqrt(ccrms/ncclip);
      disp(['Center clipped peak rms: ', num2str(ccrms)])
      range = Nstd*ccrms
      gain = (MAXSHORT-MINSHORT)/(2*range+1)
      
    case 'centercut_rms'
      idx_ccut = find(abs(normsong) > ccut);
      ccrms = std(normsong(idx_ccut));
      disp(['Center clipped rms: ', num2str(ccrms)])
      range = Nstd*ccrms
      gain = (MAXSHORT-MINSHORT)/(2*range+1)
     
    case 'rms'
      range = Nstd*songstd
      gain = (MAXSHORT-MINSHORT)/(2*range+1)
    case 'set_rms'
      gain = setrmsval/songstd
      
    case 'set_range'
      range = setrangeval;
      gain = (MAXSHORT-MINSHORT)/(2*range+1)
      
    case 'max'
      gain = .9*(MAXSHORT-MINSHORT)/(2*songmax)
      
    otherwise
      disp('Using default normalization of 1')
      gain = 1;
  end
  
  normsong = fix(sqrt(.1)*normsong*gain);
  
  % Resample if needed
  if Fsw ~= Fs
    switch resamptype
      case 'sox'
	normsongfile = [songfile, '.norm']
	resampsongfile = [songfile, '.norm.resamp']
	wavwrite(normsong/abs(MINSHORT),Fs,nbits,normsongfile)
	resampcommand = ['! sox -V -t wav ' normsongfile ' -t wav -r ' ...
	      int2str(Fsw) ' ' resampsongfile ' resample'] 
	eval(resampcommand)
	[normsong,Fsr,nbits] = wavread(resampsongfile);
	if Fsr ~= Fsw
	  disp(['Resampling failed for file: ', songfile])
	end
	normsong = normsong*abs(MINSHORT);
      case 'matlab'
	normsong = resample(normsong,Fsw,Fs,50);
      otherwise
	error('Invalid resample method.')	
    end
    normsonglen  = length(normsong)
  end
  
  pclip_idx = find(normsong > MAXSHORT);
  npclip = length(pclip_idx)
  normsong(pclip_idx) = MAXSHORT;
  nclip_idx = find(normsong < MINSHORT);
  nnclip = length(nclip_idx)
  normsong(nclip_idx) = MINSHORT;
  nclip(ifile) = npclip + nnclip
  
  disp(['Number of points clipped: ', int2str(nclip(ifile)), ', ', ...
	num2str(100*nclip(ifile)/normsonglen),'%'])
  
  % Change bit resolution if needed
  if nwbits == 32
    bitconvert = MIN18BIT/MINSHORT
    normsong = normsong*bitconvert;
  end

  % Write to file
  if songtype == 2
    normdcpsongfile = [songfile, '-10dB.norm_40k.dcp']
  else
    normdcpsongfile = [songfile, '-10dB.norm.dcp']
  end
  fidw = fopen(normdcpsongfile,'w','b')
  fprintf(fidw, 'AD_FREQ: %d Hz\n', Fsw);

  dcpsonglen = normsonglen
  maxsonglen = fix(MAXSONGDUR*Fsw)
  if dcpsonglen > maxsonglen
    dcpsonglen = maxsonglen
    disp(['Data will be truncated for file: ', normdcpsongfile])
  end
  fwrite(fidw, dcpsonglen, 'int32')

  numwrite = fwrite(fidw,normsong(1:dcpsonglen),wtype);
  if numwrite ~= dcpsonglen
    disp(['Not all of data was written for file: ', normdcpsongfile])
  end
  fclose(fidw);
  
end

disp('Done.')
