function [h_print, h_scale, h_edit, h_label] = uisongcontrols(action)

%build uicontrols
%set this up so everything is called here
%first initialize
%dont always save notefile in make-current
%make different actions (labeling/editing mutually exclusive
%fix highlight


%global variables
global h_labels
global onsets offsets labels Fs sm_win threshold min_int min_dur
global get_next get_prev SA_close interact
global soundfile path_notefile save_notes
global current_syl label_ptrdata

SA = getappdata(gcf, 'SA_Data');
h_main_amp = SA.handles.mainAmp;
h_main_labels = SA.handles.mainLabels;
do_patch = SA.params.do_patch;

if strcmp(action,'initialize')  

  %button parameters
   but_width = .03;
   but_height = .04;
   left_bord = .01;
   top_bord = .01;
   
  %control positions
   print_position = [left_bord 1-top_bord-but_height but_width but_height];
   scale_position = [left_bord 1-2*top_bord-2*but_height but_width but_height];
   edit_position = [left_bord 1-3*top_bord-3*but_height but_width but_height];
   label_position = [left_bord 1-4*top_bord-4*but_height but_width but_height];
   zoom_x_position = [left_bord 1-5*top_bord-5*but_height but_width but_height];
   zoom_y_position = [left_bord 1-6*top_bord-6*but_height but_width but_height];
   unzoom_x_position = [left_bord 1-7*top_bord-7*but_height but_width but_height];
   unzoom_y_position = [left_bord 1-8*top_bord-8*but_height but_width but_height];
   resegment_position = [left_bord 1-9*top_bord-9*but_height but_width but_height];
   move_left_position = [left_bord 1-10*top_bord-10*but_height but_width but_height];
   move_right_position = [left_bord 1-11*top_bord-11*but_height but_width but_height];
   prev_position = [left_bord 1-12*top_bord-12*but_height but_width but_height];
   next_position = [left_bord 1-13*top_bord-13*but_height but_width but_height];
   close_position = [left_bord 1-14*top_bord-14*but_height but_width but_height];
   batch_position = [left_bord 1-15*top_bord-15*but_height but_width but_height];
   psd_position = [left_bord 1-16*top_bord-16*but_height but_width but_height];

  %-------------------------------------------------------------
  %  print button
  %-------------------------------------------------------------
     h_print = uicontrol(gcf,...
                    'style','pushbutton',...
                    'units','normalized',...
                    'position', print_position,...
                    'string','p',...
		    'TooltipString', 'Print',...
                    'callback','make_current(get(0,''pointerwindow''),1); print_song');

  %-------------------------------------------------------------
  %  scale button 
  %-------------------------------------------------------------
    h_scale = uicontrol(gcf,...
                    'style','pushbutton',...
                    'units','normalized',...
                    'position',scale_position,...
                    'string','s',...
		    'TooltipString', 'Scale',...
                    'callback','make_current(get(0,''pointerwindow''),1); uiscale');

  %-------------------------------------------------------------
  % edit button 
  %-------------------------------------------------------------
    h_edit = uicontrol(gcf,...
                    'style','pushbutton',...
                    'units','normalized',...
                    'position',edit_position,...
                    'string','e',...
		    'TooltipString', 'Edit notes',...
                    'callback','uisongcontrols(''edit_note'')');

  %-------------------------------------------------------------
  % label button 
  %-------------------------------------------------------------
    h_label = uicontrol(gcf,...
                    'style','pushbutton',...
                    'units','normalized',...
                    'position',label_position,...
                    'string','l',...
		    'TooltipString', 'Label notes',...
                    'callback','uisongcontrols(''label'')');
                                  

  %-------------------------------------------------------------
  % zoom_x button 
  %-------------------------------------------------------------
    h_zx = uicontrol(gcf,...
                    'style','pushbutton',...
                    'units','normalized',...
                    'position',zoom_x_position,...
                    'string','zx',...
		    'TooltipString', 'Zoom X',...
                    'callback', 'uisongcontrols(''zoomx'')');

  %-------------------------------------------------------------
  % zoom_y button 
  %-------------------------------------------------------------
    h_zy = uicontrol(gcf,...
                    'style','pushbutton',...
                    'units','normalized',...
                    'position',zoom_y_position,...
                    'string','zy',...
		    'TooltipString', 'Zoom Y',...
                    'callback', 'uisongcontrols(''zoomy'')');

  %-------------------------------------------------------------
  % unzoom_x button 
  %-------------------------------------------------------------
    h_ux = uicontrol(gcf,...
                    'style','pushbutton',...
                    'units','normalized',...
                    'position',unzoom_x_position,...
                    'string','ux',...
		    'TooltipString', 'Unzoom X',...
                    'callback','make_current(get(0,''pointerwindow''),1); uizoom(''ux'')');

  %-------------------------------------------------------------
  % unzoom_y button          	  
  %-------------------------------------------------------------
    h_uy = uicontrol(gcf,...
                    'style','pushbutton',...
                    'units','normalized',...
                    'position',unzoom_y_position,...
                    'string','uy',...
		    'TooltipString', 'Unzoom Y',...
                    'callback','make_current(get(0,''pointerwindow''),1); uizoom(''uy'')');

  %-------------------------------------------------------------
  % resegment button 
  %-------------------------------------------------------------
    h_resegment = uicontrol(gcf,...
                    'style','pushbutton',...
                    'units','normalized',...
                    'position', resegment_position,...
                    'string','r',...
		    'TooltipString', 'Resegment notes',...
                    'callback','make_current(get(0,''pointerwindow''),2); uiresegment');

  %-------------------------------------------------------------
  % move left button 
  %-------------------------------------------------------------
    h_mv_left = uicontrol(gcf,...
                    'style','pushbutton',...
                    'units','normalized',...
                    'position', move_left_position,...
                    'string','<',...
		    'TooltipString', 'Scroll left',...
                    'callback','make_current(get(0,''pointerwindow''),1); uimove(''left'')');


  %-------------------------------------------------------------
  % move right button 
  %-------------------------------------------------------------
    h_mv_right = uicontrol(gcf,...
                    'style','pushbutton',...
                    'units','normalized',...
                    'position', move_right_position,...
                    'string','>',...
		    'TooltipString', 'Scroll right',...
                    'callback','make_current(get(0,''pointerwindow''),1); uimove(''right'')');


  %-------------------------------------------------------------
  % prev button 
  %-------------------------------------------------------------
  %make current always saves current notedata, so next saves data before opening new window
    h_prev = uicontrol(gcf,...
                    'style','pushbutton',...
                    'units','normalized',...
                    'position', prev_position,...
                    'string','<<',...
		    'TooltipString', 'Previous song file',...
                    'callback','uisongcontrols(''prev'')');


  %-------------------------------------------------------------
  % next button 
  %-------------------------------------------------------------
  %make current always saves current notedata, so next saves data before opening new window
    h_next = uicontrol(gcf,...
                    'style','pushbutton',...
                    'units','normalized',...
                    'position', next_position,...
                    'string','>>',...
		    'TooltipString', 'Next song file',...
                    'callback','uisongcontrols(''next'')');


  %-------------------------------------------------------------
  % close button 
  %-------------------------------------------------------------
  %make current always saves current notedata, so close saves data before closing window
    h_close = uicontrol(gcf,...
                    'style','pushbutton',...
                    'units','normalized',...
                    'position', close_position,...
                    'string','c',...
		    'TooltipString', 'Close',...
                    'callback','uisongcontrols(''close'')');
                    
                    
                    
  %-------------------------------------------------------------
  % batch button 
  %-------------------------------------------------------------
  
    h_batch = uicontrol(gcf,...
                    'style','pushbutton',...
                    'units','normalized',...
                    'position', batch_position,...
                    'string','b',...
		    'TooltipString', 'Batch',...
                    'callback','uisongcontrols(''batch'')');
                                        
                    
 %-------------------------------------------------------------
  % psd button 
  %-------------------------------------------------------------
  
    h_psd = uicontrol(gcf,...
                    'style','pushbutton',...
                    'units','normalized',...
                    'position', psd_position,...
                    'string','p',...                    
		    'TooltipString', 'PSD analysis',...
 		    'callback','make_current(get(0,''pointerwindow''),1); uipsdanal');

elseif strcmp(action,'label')
  
  % Get the custom pointer data.
  label_ptrdata = textptr;
  % Get the matlab version
  matver = version('-release');
  
  % make_current(get(0,'pointerwindow'),1);
  h_fig = gcbf;
  make_current(h_fig, 1);
  uistate = uisuspend(h_fig);
  % Turn off panning callbacks
  %set(h_fig, 'WindowButtonMotionFcn', '', ...
  %          'WindowButtonDownFcn', '', ...
  %          'Interruptible', 'off', ...
  %          'BusyAction', 'queue');

  % Current starting syllable is the left most one in the current view.
  nsyls = length(onsets);
  xlim = get(h_main_amp, 'XLim');
  prev_syl = nnz(xlim(1)>onsets);
  current_syl = min(prev_syl+1, nsyls)
  
  subplot(h_main_amp)
  %[h_amp_plot, h_labels]=disp_song;
  title('LABEL MODE! (q to quit)','color',[1,1,0]);      
  set(h_labels(current_syl),'color',[1 0 0]);

  % KeyPressFcn callback should remove these callbacks and reset
  % default ones when 'q' is pressed.
  %set(0, 'CurrentFigure', h_fig);

  if matver <= 12.1
    set(h_fig,'WindowButtonMotionFcn', 'label_select_clbk_m61(''label_pointer'')', ...
              'WindowButtonDownFcn', 'label_select_clbk_m61(''select'')');
    
    set(h_fig,'KeyPressFcn', 'label_keypress_clbk_m61(uistate)', ... 
              'Interruptible', 'off', 'BusyAction', 'queue');
  else
    set(h_fig,'WindowButtonMotionFcn', {@label_select_clbk, 'label_pointer', onsets, offsets}, ...
              'WindowButtonDownFcn',{@label_select_clbk, 'select', onsets, offsets});
    
    set(h_fig,'KeyPressFcn', {@label_keypress_clbk, onsets, offsets, uistate}, ... 
              'Interruptible', 'off', 'BusyAction', 'queue');
  end
  %set(h_main_amp,'ButtonDownFcn', {@label_select_clbk, onsets, offsets});
  % This doesn't work because the patch is in the way!
  %set(h_main_labels,'ButtonDownFcn', {@label_select_clbk, onsets, offsets});

elseif strcmp(action,'prev')
  make_current(get(0,'pointerwindow'),1);
   
  %save data from the current window if save flag is set to 1
  if save_notes==1
    if (ischar(soundfile) & ~isempty(soundfile) & ~isspace(soundfile) )
      note_file=[soundfile,'.not.mat'];
      save_data(note_file, path_notefile, Fs, onsets, offsets, labels, threshold, min_int, min_dur, sm_win); 
    end
  end
  %set prev flag
  get_prev = 1; 

elseif strcmp(action,'next')
  make_current(get(0,'pointerwindow'),1);
    
  %save data from the current window if save flag is set to 1
  if save_notes==1
    if (ischar(soundfile) & ~isempty(soundfile) & ~isspace(soundfile) )
      note_file=[soundfile,'.not.mat'];
      save_data(note_file, path_notefile, Fs, onsets, offsets, labels, threshold, min_int, min_dur, sm_win); 
    end
  end
  %set next flag
  get_next = 1; 

elseif strcmp(action,'close')
   make_current(get(0,'pointerwindow'),1);
   %save data from the current window if save_notes=1
   if save_notes==1
     if (ischar(soundfile) & ~isempty(soundfile) & ~isspace(soundfile) )
       note_file=[soundfile,'.not.mat'];
       save_data(note_file, path_notefile, Fs, onsets, offsets, labels, threshold, min_int, min_dur, sm_win); 
     end
   end
   SA_close = 1;
   %close window
   %close(gcf);
   
elseif strcmp(action,'batch')
   make_current(get(0,'pointerwindow'),1);
   %save data from the current window if save_notes=1
   if save_notes==1
     if (ischar(soundfile) & ~isempty(soundfile) & ~isspace(soundfile) )
       note_file=[soundfile,'.not.mat'];
       save_data(note_file, path_notefile, Fs, onsets, offsets, labels, threshold, min_int, min_dur, sm_win); 
     end
   end
   interact = batch_setup;

elseif strcmp(action,'zoomx')
  make_current(get(0,'pointerwindow'),1); 
  uizoom('x', do_patch);

elseif strcmp(action,'zoomy')
  make_current(get(0,'pointerwindow'),1); 
  uizoom('y', do_patch);

elseif strcmp(action,'edit_note')
  make_current(get(0,'pointerwindow'),1); 
  edit_note(do_patch);

end

