function [left,right,button,flag]=get_xrange(varargin)

%allows user to select range of x axis values with mouse from current axis
%select left value with left button
%select right value with middle or right button
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
ymin=ylim(1);
ymax=ylim(2);

%default values
left=xmin;
right=xmax;
button=1;
flag=0;  % no change state

%plot default lines

h_left=line([left left],[ymin ymax],'color',[0,1,0],'LineStyle','--');
h_right=line([right right],[ymin ymax],'color',[1,0,0],'LineStyle','--');

if do_patch
  h_patch = patch('xdata',[left left right right], ...
                  'ydata',[ymin ymax ymax ymin], ...
                  'FaceAlpha', 0.2, 'FaceColor', 'y', ...
                  'LineStyle', 'none');
end

while ~isempty(button) & ((button == 1) | (button == 2) | (button == 3))

  [xval,y,button]=myginput(1,h_axis)
  
  if ~isempty(button)   
    if (button == 1)
      left=max(xval,xmin);
      left=min(left,xmax);
      set(h_left,'xdata',[left left]) 
      if do_patch
        if left < right
          set(h_patch,'xdata',[left left right right], ...
                      'ydata',[ymin ymax ymax ymin], 'Visible','on');
        else
          set(h_patch,'Visible','off')
        end
      end
    elseif (button == 2) | (button == 3)
      right=min(xval,xmax);
      right=max(right,xmin);
      set(h_right,'xdata',[right right])     %move line
      if do_patch
        if left < right
          set(h_patch,'xdata',[left left right right],'ydata',[ymin ...
                              ymax ymax ymin],'Visible','on');
        else
          set(h_patch,'Visible','off')
        end
      end
    end
  end

end
 

%if quiting, or illegal values, or no change leave values unchanged
%and set flag to 0

if (left == xmin) & (right == xmax)
   flag  = 0;
elseif (left >= right) | strcmp(char(button),'q')   
   left=xmin
   right=xmax
   flag = 0;
else
   flag = 1
end

%delete lines

delete(h_left)
delete(h_right)
if do_patch
  delete(h_patch)
end


     
