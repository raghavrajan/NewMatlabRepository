function uisonganal_copy(read_flag)

%optional argument read_flag: if read_flag=='ro' (read only) then don't try to write data 
   if nargin == 0
     read_flag='rw';   %read and write by default
   end   

%some global variables   

   %values for display
   global d_song time_d_song     %values used by disp_song
   global spect time_spect freq_spect
   global threshold smooth Fs
   global h_main h_main_spect h_main_amp
   global win_percent
   global pr_left_margin pr_bottom_margin
   global main_win_code
   global Fs
   global spect_floor spect_ceil spect_range
   global soundfile
   
   
   %control flow
   global get_next interact batch_disp save_notes save_spect save_filt batch_print

   %values used in segmentation
   global threshold min_int min_dur sm_win
   global onsets offsets labels
   global smooth filtsong

   %file data and pathnames
   global soundfile path_notefile path_filtfile path_spectfile

   %indexes into main figure userdata
   %these just facilitate mneumonics of storage and retrieval of values in userdata
   %values are set below and used by programs that retrieve from userdata
   global imain_win_code idata_status ih_main_amp ih_main_spect ih_fname 
   global ih_path_songfile ih_path_notefile ih_path_filtfile ih_path_spectfile ih_amp_plot  
   
   
%some default values

   %identification code for figure containing spectrogram and amplitude plots
   % stored as first value of figure userdata
   main_win_code = -101;

   %indexes for userdata storage and retrieval
   %code identifying window type (number)
     imain_win_code = 1;
   %code indicating status of data for current figure
   %see make_cuurent for meaning of status codes
     idata_status = 2;
   %handles to amp and spect axes
     ih_main_amp = 3;
     ih_main_spect = 4;
   %handles to hidden text containing filename and path strings
     ih_fname = 5;
     ih_path_songfile = 6;
     ih_path_notefile = 7;
     ih_path_filtfile = 8;
     ih_path_spectfile = 9;
   %handles to amplitude plot
     ih_amp_plot = 10;
     
     
   %interactive or batch mode?
   interact = 1;
   
   %display data during batch processing?
   batch_disp = 1;
   
   %print in batch mode?
   batch_print = 0;
   
   %save notefiles, filtered song, spectrogram?
   if strcmp(read_flag,'ro')
     save_notes = 0;
     save_filt = 0;
     save_spect = 0;
   elseif strcmp(read_flag,'rw')   
     save_notes = 1;
     save_filt = 1;
     save_spect = 1;
   end

   
   get_next = 0;

   %main window handle, starts unset
   h_main = [];
   
   %for bandpass
   F_low=300;
   F_high=8000;

   %window length for smoothing prior to segmentation (in ms)
   sm_win=2;
 
   %threshold for segmentation, a value of 0 causes threshold to be calculated
   threshold = 0;     %is there a current threshold? 
      
   %d_Fs is resample rate for more rapid display of song
   d_Fs=1000;
   
   %criteria for intervals and notes (minimum durations in ms)
   min_int=5;
   min_dur=10;
   
   %values for calculating spectrogram
   spect_win_dur=16;
   spect_overlap = .50;  %percentage of overlap of specgram window
   
   
   %initial floor, ceiling and contrast for display of spectrogram
   spect_floor = -55;  %floor cutoff in dB re max
   spect_ceil = 0;   %ceil saturation in db re max
   spect_range = 2.0; % percent of range of grayscale used to display data between floor and ceil
   
      
   %win_percent sets the amount by which window is shifted left or right 
   win_percent = 0.80;
   
   %default values for printing
   pr_orient = 'landscape';  %paper oreintation
   pr_papertype = 'usletter'; %papertype 
   pr_y_paperinches = 6; %number of inches for height of printed figure
   pr_x_per_paperinch = 600; %number of units(ms) per inch of printed figure
   pr_left_margin = .0;
   pr_bottom_margin = .0;
            

%get name of metafile containing songfile names

meta_fid = -1;
metafile = 0;
   disp('select batchfile');
   [metafile, pathname]=uigetfile('*batch*','select batchfile')
   meta_fid=fopen([pathname, metafile]);
   if meta_fid == -1
      disp('cannot open file' )
      disp (metafile)
   end

%set paths for storage and retrieval, for now default to pathname
path_songfile = pathname;
path_notefile = pathname;
path_filtfile = pathname;
path_spectfile = pathname;

%get filetype, later versions of read program can be smarter about this

disp('What is type of sound file? [w]')
disp(' b = binary from mac')
disp(' w = wavefile (i.e. cbs)')
disp(' d = dcpfile')
disp(' f = foosong/gogo')
disp(' o = observer file (last/song channel)')
disp(' o1r = observer file (second to last)')


filetype = 'null';
while strcmp(filetype,'null')
  temp=input(' ','s');
  if (strcmp(temp,'b') | strcmp(temp,'B') | strcmp(temp,'bin') | strcmp (temp, 'binary'));
     filetype = 'b';
  elseif strcmp(temp,'w')  | strcmp(temp,'wav') | strcmp(temp,'W') | isempty(temp)
     filetype = 'w';
  elseif strcmp(temp,'d')  | strcmp(temp,'dcp') | strcmp(temp,'D')
     filetype = 'd';  
  elseif strcmp(temp,'f')  | strcmp(temp,'F')
     filetype = 'f';  
  elseif strcmp(temp,'o')  | strcmp(temp,'O')
     filetype = 'obs0r';
  elseif strcmp(temp,'o1r')  | strcmp(temp,'O1r')
     filetype = 'obs1r';  else
     disp('Unacceptable! Pick again')
  end
end

%get default sample rate

default_Fs=input('What is sample rate (samples/sec)? [32000] ');
if isempty(default_Fs); default_Fs = 32000; end

%main program: cycle through sound files until end or quit command

while 1
   %reset values of onsets, offsets, labels, filtered song, spectrogram
   onsets = [];
   offsets = [];
   labels = [];
   filtsong = [];
   spect = [];
      
   %values for threshold, Fs, min_int, min_dur, sm_win will be maintained
   %unless new values are read in from a note file
     
   %get soundfile name
   soundfile = fscanf(meta_fid,'%s',1)
   if isempty(soundfile)
      disp('End of songfiles')
      break
   end
%   [path_songfile,soundfile,ext,ver] = fileparts(soundfile)   

   %if soundfile name ends in '.filt' then strip '.filt' from the end
   % so that batch files which list filtfiles can be read
   filt_idx=findstr('.filt',soundfile);
   if ~isempty(filt_idx)
     soundfile=soundfile(1:filt_idx(length(filt_idx))-1);
   end
   
   %get note_file name
    note_file=[soundfile,'.not.mat'];  
   %if notefile exists, get it
  if exist([path_notefile, note_file])
     disp('loading notefile...');     
     load([path_notefile, note_file]);
     default_Fs = Fs;   %if  Fs was read, it is the new default value
  end

   %if file was loaded then: 
   % Fs, threshold, min_int, min_dur, sm_win, onsets, offsets and labels should be defined by notefile values
   % otherwise onsets, offsets and labels are empty, but other values are set at defaults
           
   
   %if filtered song exists, get it, otherwise read and filter song, save filtered song if flag is set
   %note: if filtsong exists, it is unneccessary to have actual songfile around    
   filt_file=[path_filtfile, soundfile, '.filt'];
   if exist(filt_file)     
     disp('loading filtered song...');
     [filtsong, Fs] = read_filt(filt_file);
     default_Fs = Fs;   %if Fs was read, it is the new default value (takes precedence over notefile if discrepant
   else
   
     [rawsong,Fs]=soundin_copy(path_songfile, soundfile, filetype);
     disp('filtering song...');
     %unless Fs has been read from the file use the default value which was
     %set by the user or previously read in from a note file
     if Fs == -1; Fs = default_Fs; end
     %filter soundfile
     filtsong=bandpass(rawsong,Fs,F_low,F_high);
     %save soundfile if flag is set
     if save_filt == 1     
        disp('saving filtered song...')
        write_filt(filt_file, filtsong, Fs);
     end
   end

   %if spectrogram is going to be saved or displayed
   % then see if it already exists, if not, calculate it and save if flag is set
   if (batch_disp == 1) | (save_spect == 1);
      spect_file=[path_spectfile, soundfile,'.spect'];
      if exist(spect_file)     
        disp('reading spectrogram...')
        [idx_spect, nfft, spect_win, noverlap, t_min, t_max, f_min, f_max] = read_spect(spect_file);
        time_spect = [t_min, t_max];
        freq_spect = [f_min, f_max];     
      else
        %calculate spectrogram
        disp('calculating spectrogram...');
        %first calculate nfft and noverlap
        nfft=round(Fs*spect_win_dur/1000);
	nfft = 2^nextpow2(nfft);
        spect_win = hanning(nfft);
        noverlap = round(spect_overlap*length(spect_win)); %number of overlapping points       
        %now calculate spectrogram
        [spect, freq, time] = specgram(filtsong, nfft, Fs, spect_win, noverlap);
        idx_spect=scale_spect(spect);  %calculate index array for spectrogram
        f_min = freq(1);
        f_max = freq(length(freq));
        freq_spect = [f_min, f_max];
        t_min = time(1)*1000; %convert to ms
        t_max = time(length(time))*1000; %convert to ms
        %adjust time axis for spectrogram offset (1/2 window duration in ms)
        t_min = t_min + 0.5*nfft*1000/Fs;  
        t_max = t_max + 0.5*nfft*1000/Fs;  
        time_spect = [t_min, t_max];                
        %save spectrogram if flag is set        
        if save_spect == 1     
          disp('saving spectrogram...')
          write_spect(spect_file, idx_spect, nfft, spect_win, noverlap, t_min, t_max, f_min, f_max);
        end
      end      
   end   

   disp('calculating power and smoothing...')
   %calculate square of signal: proportional to power
   squared_song = filtsong.^2;
   
   %smooth the rectified song
   len=round(Fs*sm_win/1000);                      
   h=ones(1,len)/len;
   smooth=conv(h, squared_song);
   offset=round((length(smooth)-length(filtsong))/2); %get rid of convolution induced offset
   smooth=smooth(1+offset:length(filtsong)+offset);
   
   %if threshold is undefined, calculate threshold and segment: will only be true for first file
   if threshold == 0
      disp('calculating default threshold...')
      %get threshold value
      len1=round(Fs*10/1000);
      [m, i]=max(smooth(len1:length(smooth)-len1));
      peak=sum(smooth(i-len1:i+len1))/(2*len1+1);
      threshold=round(peak)/1000;
   end   
   
   
   %if values were not loaded for onsets offsets and labels, calculate them, and save them if flag is set
   if isempty(onsets)   
      disp('segmenting song...')
      [onsets, offsets] = segment(smooth, Fs, min_int, min_dur, threshold);
      %set label values
      labels= num2str(zeros(size(onsets)))';
      %save note data if flag is set
      if save_notes == 1; 
        save_data(note_file, path_notefile, Fs, onsets, offsets, labels, threshold, min_int, min_dur, sm_win);
      end 
   end  
   
   
   %fix for matlab 5.x compatability
   
  labels=makerow(labels);
  onsets=makecol(onsets);
  offsets=makecol(offsets);      

 %display only if flag is set
 if batch_disp == 1 | batch_print == 1
   
   %resample song for more rapid display
   disp('resampling for display...');
   step=round(Fs/d_Fs);
   d_Fs=Fs/step;
   d_song=smooth(1:step:length(smooth));

   
   %set initial values of display variables  
   time_d_song=[0:length(d_song)-1]*1000/d_Fs;  %vector for time axis
   g_xmin=time_d_song(1);                                    %initially display everything
   g_xmax=time_d_song(length(time_d_song));
   g_ymin=0;
   g_ymax=max(d_song)*1.2;
   % The following wierdness is a workaround because of matlab bug for image display
   %see disp_spect for a bit more info
   %basic idea is to invert the spectrogram before displaying and to label y axis with negative 
   %frequency values
   g_fmin=(-1)*F_high;
   g_fmax=0;  
     
   %display song
   disp('displaying song...')
   %could retain old figure here
   if isempty(h_main);
      h_main=figure;
   end 
   %see if this helps w/ gradual slowing of program as files are processed
   %works but doesn't preserve window size
   %delete(h_main);
   %h_main=figure;
   
   h_main_amp=subplot(2,1,2);
   %store the global values of data limits with axes
   set(h_main_amp,'userdata',[g_xmin g_xmax g_ymin g_ymax]);
   %set the current display to show all the data
   set(h_main_amp,'xlim',[g_xmin g_xmax]);
   set(h_main_amp,'ylim',[g_ymin g_ymax]);
   h_amp_plot = disp_song;
   
   %display spectrogram
   disp('displaying spectrogram...')
   h_main_spect=subplot(2,1,1);
   set(h_main_spect,'ydir', 'reverse');
    %initially display all data
   set(h_main_spect,'xlim',[g_xmin g_xmax]);
   set(h_main_spect,'ylim',[g_fmin g_fmax]);
   %store max and min data values in user data
   set(h_main_spect,'userdata',[g_xmin g_xmax g_fmin g_fmax]);
   title(soundfile);
   disp_idx_spect(idx_spect, time_spect, freq_spect, spect_floor, spect_ceil, spect_range);
   subplot(h_main_amp);

   %set some printing defaults
   set_print_vals(h_main, pr_orient, pr_papertype, pr_x_per_paperinch, g_ymax-g_ymin); %y value is superceded on next line
   temp = get(h_main, 'paperposition');
   temp(4) = pr_y_paperinches;
   set(h_main, 'paperposition', temp);
      
   %reset data status for any open windows
   windows = get(0,'children');
   for i = 1:length(windows)
        userdata = get(windows(i),'userdata');
        if length(userdata > 0)
             if userdata(imain_win_code) == main_win_code;
                 userdata(idata_status) = 0;
                 set(windows(i),'userdata',userdata);
             end
        end
   end

   %store file and pathnames as hidden text
   %store these with spectrogram axes for now: note this will be a problem if spectrogram axis
   %is cleared
   subplot(h_main_spect);
   h_soundfile = text(-100, -100, soundfile, 'visible', 'off');
   h_path_songfile = text(-100, -100, path_songfile, 'visible', 'off');                       
   h_path_notefile = text(-100, -100, path_notefile, 'visible', 'off');   
   h_path_filtfile = text(-100, -100, path_filtfile, 'visible', 'off'); 
   h_path_spectfile = text(-100, -100, path_spectfile, 'visible', 'off'); 
    
   %set userdata values
   userdata(imain_win_code) = main_win_code;
   userdata(idata_status) = 3;  %see make_current for meaning of code
   userdata(ih_main_amp) = h_main_amp;
   userdata(ih_main_spect) = h_main_spect;
   userdata(ih_fname) = h_soundfile; 
   userdata(ih_path_songfile) = h_path_songfile;  
   userdata(ih_path_notefile) = h_path_notefile;
   userdata(ih_path_filtfile) = h_path_filtfile;
   userdata(ih_path_spectfile) = h_path_spectfile;
   userdata(ih_amp_plot) = h_amp_plot;    
   
   subplot(h_main_amp);   
   
   %store info in userdata for retrieval of handles etc
   set(h_main,'userdata',userdata);
  
   %build ui controls: sets up all the interface controls in current figure
   uisongcontrols('initialize');
   uispectcontrols('initialize');  
   
   %if in batch mode and want to display, need to pause here
   if interact == 0 & batch_disp == 1
      pause(.3);
   end
  
  end   %end of building window if in display mode
  
  %if batch_print is set, print window
  if batch_print == 1
     %set size of current window to fit page
       set(gcf,'paperorientation', pr_orient);
     %save unit values and set to inches
       t_paperunits = get(gcf,'PaperUnits');
       set(gcf,'PaperUnits','inches');
       t_fig_units = get(gcf,'Units');
       set(gcf,'Units','inches');set(gcf,'units','inches');
     %set figure on screen to match page size
       position = get(gcf,'position');
       position(3)=11;
       position(4)= pr_y_paperinches;
       set(gcf,'position',position);
     %set horizontal scale of current window to desired value
       set_scale_vals(gcf, pr_x_per_paperinch, pr_y_paperinches);
     %set paperposition values
       position(1) = pr_left_margin;
       position(2) = pr_bottom_margin;
       set(gcf,'paperposition',position);
     %print & shift until done with spectrogram
       more = 1;
       while more == 1
         print -dps2 -noui
         subplot(h_main_amp);
         more = move_right(win_percent);
         subplot(h_main_spect);
         move_right(win_percent);
         subplot(h_main_amp);
       end
     %reset units  
      set(gcf,'PaperUnits',t_paperunits);
      set(gcf,'Units',t_fig_units);
  end
  
  
  %sit idle until need to build another window or in batch mode
  while get_next==0 & interact == 1
     pause(.5)
  end
  
  %reset
  get_next = 0;
  
   
end

fclose(meta_fid);
