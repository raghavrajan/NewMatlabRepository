function autocorr_ff(batch_in,note)

%function psdanal_ff pulls out all examples of the sound defined by string 'note from all of the 
%songs in the file  of notefiles 'batch_in' and calculates the ff of a 16ms segment of the
%syllable, which is determined by the user
%the user enters the start of the segment by entering either the ms_from_start of the syllable
%or the % in the syllable

%this function uses parts of pick_notes and get_songs to pull out the relevant 16ms segment of the
%syllable of interest, calculates the psd of that segment, and then calculates the auto-covariance
%it looks for the distance, in frequency, between the 0th lag and the first peak

%this function assumes that the file type is 'filt' b/c it requires that notes are already labeled 
%using uisonganal

%this function assumes that you are in the subdirectory w/ the data files and the .filt files




sample_size=16;
k=1;
freq=cell(100,2);
savedata=0;         %do not save the data unless specified
    
    
%make sure note is a string
note=char(note);


%align by the start of the note if align_start==1
align_start=1;


%open batch_file
meta_fid=fopen([batch_in]);
if meta_fid==-1|batch_in==0
        disp('cannot open file')
        disp(batch_in)
end

%what part of the syllable to analyze?
syl_segment=input('1) percent from start or 2) ms from start?  ');

if syl_segment==1
    percent_from_start=input('percent from start?  ');
else if syl_segment==2
    ms_from_start=input('ms from start?   ');
end    

savedata=input('Do you want to save the values calculated?  1)yes   2) no ');
if savedata==1
    temp=input('Save as what?   ','s')
    jj=['',temp,'_ff'];
end




while 1
       %get notefile name
       notefile=fscanf(meta_fid,'%s',1);
       %end when there are no more notefiles
       if isempty(notefile);
           break
       end
       
       %if notefile exists, get it
       if exist([notefile])
           load(notefile);          %Fs, onsets, offsets and labels are defined
       else disp(['cannot find ',notefile])
       end
       
       %count the number of matches in this file
       matches=findstr(note,labels);   %returns the index for the note in the labels vector
       num_matches=length(matches);
       
       %get soundfile name
       end_file=findstr('.not.mat',notefile);
       rootfile=notefile(1:end_file-1);
       soundfile=[rootfile,'.filt'];
       
       %get the sound
       for i=1:num_matches  %for all occurrencs of the syllable
            labels=makerow(labels);
            start_time=onsets(matches(i));
            start_time=start_time/1000;         %convert onsets and offsets to seconds
            end_time=offsets(matches(i));
            end_time=end_time/1000;
            
            %get rawdata
            [rawsong,Fs]=soundin('',soundfile,'filt');
            
            %get the desired note
            sound=rawsong(start_time*32000:end_time*32000);
            
            %take a 16ms slice through the note
            if syl_segment==1
                note_length=end_time-start_time;        %note_length is in seconds
                seg_start_time=start_time+(note_length*(percent_from_start/100));
                seg_end_time=seg_start_time+.016;
            else if syl_segment==2                       %(ms from the start)
                seg_start_time=start_time+(ms_from_start/1000);     %(in seconds)
                seg_end_time=(seg_start_time+.016);                 %(16 ms segment) 
            end  
            
            newstarttime=seg_start_time*32000;
            newendtime=seg_end_time*32000;
            note_segment=rawsong(newstarttime:newendtime);
            
            %analyze the psd of selected song segment
            %nfft=8192;
            %window=8192;
            %[Pxx,freq]=psd(note_segment,nfft,Fs,window);
            %amp=sqrt(Pxx);
            
            
            %calculate the auto-covariance of the psd of the song segment
            autocorr=xcov(note_segment);
            [maxvalue,loc]=max(autocorr(length(note_segment+3):length(note_segment+100));
            loc=loc+3;
            ff=Fs/loc
            
            %save values in a cell array w/ 2 elements:  songname and fundfreq
            ffreq{k,1}=soundfile;
            ffreq{k,2}=ff;
            k=k+1;
         
            
            
            if savedata==1
            save (jj,'ffreq')         %save the names and ff in a mat file called 'ff_data'
            end
            
        end     
            
           
    end %for loop

end     %while loop

fclose(meta_fid);

end
