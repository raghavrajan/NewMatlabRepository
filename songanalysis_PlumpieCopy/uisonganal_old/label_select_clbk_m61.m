function label_select_clbk_m61(action, varargin)

global h_main h_main_amp h_main_labels h_labels
global current_syl label_ptrdata onsets offsets

% disp(action);

switch action
  case 'label_pointer'
    % disp('Got to label_pointer');
    % oldptr = get(gcbf,'Pointer');
    POINTER_IN_AXES = 0;
    fig = get(0,'PointerWindow'); 
    % Look for quick exit
    if fig==0,
      return
    end
    
    saveunits = get(gcbf,'Units');
    set(gcbf,'Units','pixels');
    figpos = get(gcbf,'Position');
    ha = findobj(gcbf,'Type','axes','Tag','SAAxis');
    for i=1:length(ha)
      ptrloc= get(0,'PointerLocation');
      ptrloc = (ptrloc - figpos(1,1:2))./figpos(1,3:4);
      % If pointer is in figure then check if it is in an axis
      if (ptrloc(1) <= 1 & ptrloc(1) >= 0 & ptrloc(2) <= 1 & ptrloc(2) >= 0)
	ax = ha(i);
	saveaxesunits = get(ax,'Units');
	set(ax,'Units','normalized');
	axispos = get(ax,'Position');
	if (axispos(1) <= ptrloc(1) & ptrloc(1) <= axispos(1)+axispos(3) & ...
	      axispos(2) <= ptrloc(2) & ptrloc(2) <= axispos(2)+axispos(4))	  
	  % disp('Got to pointer change');
	  POINTER_IN_AXES = 1;
	  set(ax,'Units',saveaxesunits);
	  break;
	else
	  set(ax,'Units',saveaxesunits);
	end % if	
      end
    end

    if POINTER_IN_AXES
      set(gcbf,'Pointer','custom','PointerShapeCData', label_ptrdata, 'PointerShapeHotSpot',[1 1])
    else
      set(gcbf,'Pointer','arrow');
    end
    set(gcbf,'Units',saveunits);
    
  case 'select'
    POINTER_IN_AXES = 0;
    fig = get(0,'PointerWindow'); 
    % Look for quick exit
    if fig==0,
      return
    end    
    
    saveunits = get(gcbf,'Units');
    set(gcbf,'Units','pixels');
    figpos = get(gcbf,'Position');
    ha = findobj(gcbf,'Type','axes','Tag','SAAxis');
    for i=1:length(ha)
      ptrloc= get(0,'PointerLocation');
      ptrloc = (ptrloc - figpos(1,1:2))./figpos(1,3:4);
      % If pointer is in figure then check if it is in an axis
      if (ptrloc(1) <= 1 & ptrloc(1) >= 0 & ptrloc(2) <= 1 & ptrloc(2) >= 0)
	ax = ha(i);
	saveaxesunits = get(ax,'Units');
	set(ax,'Units','normalized');
	axispos = get(ax,'Position');
	if (axispos(1) <= ptrloc(1) & ptrloc(1) <= axispos(1)+axispos(3) & ...
	      axispos(2) <= ptrloc(2) & ptrloc(2) <= axispos(2)+axispos(4))	  
	  POINTER_IN_AXES = 1;
	  set(ax,'Units',saveaxesunits);
	  break;
	else
	  set(ax,'Units',saveaxesunits);
	end % if	
      end
    end

    if POINTER_IN_AXES
      pt = get_currentpoint(ax);
      new_syl = max(nnz(pt(1,1)>onsets),1);
      if new_syl ~= current_syl
        % Unhighlight current syllable
        set(h_labels(current_syl),'color',[0 0 0]);         
        current_syl = new_syl;
        
        ha = findobj(h_main,'Type','axes','Tag','SAAxis');
        xlim = get_xlim(h_main_amp); 
        % If running off the edge of the display, center display on the
        % current selection.
        if onsets(current_syl) < xlim(1) | offsets(current_syl) > xlim(2)
          xrange = xlim(2) - xlim(1);
          xshift = onsets(current_syl) - (xlim(1) + xrange/2);
          xlim_new = xlim + xshift;
          for i=1:length(ha)
            set(ha(i),'XLim',xlim_new);
          end      
        end
        
        % Highlight the current syllable
        set(h_labels(current_syl),'color',[1 0 0]);         
        drawnow
      end
    end
    set(gcbf,'Units',saveunits);

end


function p = get_currentpoint(ax)
%GET_CURRENTPOINT Return equivalent linear scale current point
p = get(ax,'currentpoint'); p = p(1,1:2);
if strcmp(get(ax,'XScale'),'log'),
   p(1) = log10(p(1));
end
if strcmp(get(ax,'YScale'),'log'),
   p(2) = log10(p(2));
end

function xlim = get_xlim(ax)
%GET_XLIM Return equivalent linear scale xlim
xlim = get(ax,'xlim');
if strcmp(get(ax,'XScale'),'log'),
   xlim = log10(xlim);
end

function ylim = get_ylim(ax)
%GET_YLIM Return equivalent linear scale ylim
ylim = get(ax,'ylim');
if strcmp(get(ax,'YScale'),'log'),
   ylim = log10(ylim);
end


