function [bottom,top,button,flag]=get_yrange(varargin)

%allows user to select range of y axis values with mouse using current axis
%select top value with left button
%select bottom value with middle or right button
%hit any key (except q) when done
%q to quit without changing values
%the key value is returned as button
%flag is set to 1 if legal changes were made

if nargin == 0
  do_patch = 1;
  h_axis = gca;
elseif nargin == 1
  do_patch = varargin{1};
  h_axis = gca;
elseif nargin > 1
  do_patch = varargin{1};
  h_axis = varargin{2};
end

%initial values from current axes
xlim=get(gca,'xlim');
ylim=get(gca,'ylim');
xmin=xlim(1);
xmax=xlim(2);
ymin=ylim(1)
ymax=ylim(2)

%default values
top=ymax;
bottom=ymin;
button=1;
flag=0; %no change state

%plot default lines

h_top=line([xmin xmax],[top top],'color',[0,1,0],'LineStyle','--');
h_bottom=line([xmin xmax],[bottom bottom],'color',[1,0,0],'LineStyle','--');

if do_patch
  h_patch = patch('xdata',[xmin xmin xmax xmax],'ydata',[bottom top ...
                      top bottom], 'FaceAlpha', 0.2, 'FaceColor', 'y', 'LineStyle', 'none');
end

while ~isempty(button) & ((button == 1) | (button == 2) | (button == 3))

  [x,yval,button]=myginput(1,h_axis)
  
  if ~isempty(button)
    if (button == 1)
      top=max(yval,ymin);
      top=min(top,ymax);
      set(h_top,'ydata',[top top]);
      if do_patch
        if top > bottom
          set(h_patch,'xdata',[xmin xmin xmax xmax],'ydata',[bottom top ...
                              top bottom],'Visible','on');
        else
          set(h_patch,'Visible','off');
        end
      end
    elseif (button == 2) | (button == 3)
      bottom=min(yval,ymax);
      bottom=max(bottom,ymin);
      set(h_bottom,'ydata',[bottom bottom]);     %move line
      if do_patch
        if top > bottom
          set(h_patch,'xdata',[xmin xmin xmax xmax],'ydata',[bottom top ...
                              top bottom],'Visible','on');
        else
          set(h_patch,'Visible','off');
        end
      end
    end
  end

end

%if quiting, or illegal values, or no change leave values unchanged
%and set flag to 0
if (bottom == ymin) & (top == ymax)
   flag  = 0;
elseif (bottom >= top) | strcmp(char(button),'q')
   bottom=ymin;
   top=ymax;
   flag=0;
else
   flag=1;
end

%delete lines

delete(h_top)
delete(h_bottom)
if do_patch
  delete(h_patch)
end
     
