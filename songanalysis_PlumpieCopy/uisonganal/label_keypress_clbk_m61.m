function label_keypress_clbk_m61(uistate, varargin)

global h_main h_main_amp h_main_labels h_labels
global current_syl labels onsets offsets

key = double(get(h_main,'CurrentCharacter'));
ch = char(key);

if isempty(key) | isempty(ch) | ~ischar(ch)
  disp('Key empty or not a char.');
  return;
end

if strcmp(ch,'q')
  % Unhighlight current syllable
  set(h_labels(current_syl),'color',[0 0 0]);         
  subplot(h_main_amp);
  title('');
  % remove callbacks
  set(h_main_amp,'ButtonDownFcn', '');
  set(h_main_labels,'ButtonDownFcn', '');
  set(h_main,'KeyPressFcn', '');
  % Restore default callbacks
  uirestore(uistate);  
  %set(h_main,'WindowButtonMotionFcn', ...
  %           'centercb(''center_pointer'')', ...
  %           'WindowButtonDownFcn','centercb(''center'')', ...
  %           'Interruptible', 'on');
  return;
end

nsyls = length(onsets);

if key == 29 | key == 13
  %disp('Next syllable')
  % go to next syllable; 
  % move display to the right if syllable is not visible
  
  % Unhighlight current syllable
  set(h_labels(current_syl),'color',[0 0 0]);         
  current_syl = min(current_syl+1,nsyls);
  
  xlim = get(h_main_amp,'xlim'); 
  if offsets(current_syl) > xlim(2)      
    % Center on start of next syllable
    ha = findobj(h_main,'Type','axes','Tag','SAAxis');
    xrange = xlim(2) - xlim(1);
    xshift = onsets(current_syl) - (xlim(1) + xrange/2);
    xlim_new = xlim + xshift;
    for i=1:length(ha)
      set(ha(i),'XLim',xlim_new);
    end      
  end    
  
  % Highlight the current syllable
  set(h_labels(current_syl),'color',[1 0 0]);         
  drawnow;
elseif key == 28 | key == 8 | key == 127
  %disp('Previous syllable')
  % go to previous syllable
  % move display to the left if syllable is not visible
  
  % Unhighlight current syllable
  set(h_labels(current_syl),'color',[0 0 0]);         
  current_syl = max(current_syl-1,1);
  
  xlim = get(h_main_amp,'xlim'); 
  if onsets(current_syl) < xlim(1)      
    % Center on start of previous syllable
    ha = findobj(h_main,'Type','axes','Tag','SAAxis');
    xrange = xlim(2) - xlim(1);
    xshift = onsets(current_syl) - (xlim(1) + xrange/2);
    xlim_new = xlim + xshift;
    for i=1:length(ha)
      set(ha(i),'XLim',xlim_new);
    end      
  end
  
  % Highlight the new current syllable
  set(h_labels(current_syl),'color',[1 0 0]);         
  drawnow;
else
  % label current syllable with current character.
  set(h_labels(current_syl),'String', ch);
  labels(current_syl)=ch;
  
  % Jump to the next syllable and make it current
  %disp('Next syllable');
  
  % Unhighlight current syllable
  set(h_labels(current_syl),'color',[0 0 0]);         
  current_syl = min(current_syl+1,nsyls);
  
  xlim = get(h_main_amp,'xlim'); 
  if offsets(current_syl) > xlim(2)      
    % Center on start of next syllable
    ha = findobj(h_main,'Type','axes','Tag','SAAxis');
    xrange = xlim(2) - xlim(1);
    xshift = onsets(current_syl) - (xlim(1) + xrange/2);
    xlim_new = xlim + xshift;
    for i=1:length(ha)
      set(ha(i),'XLim',xlim_new);
    end      
  end    
  
  % Highlight the current syllable
  set(h_labels(current_syl),'color',[1 0 0]);         
  drawnow;
end


  