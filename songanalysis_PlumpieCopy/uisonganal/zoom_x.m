function [left, right, button, flag] = zoom_x(varargin)

%zoom in to x axes of current axes
%select left and right boundaries with left and right mouse buttons
%hit any key when done (q to quit without updating)
%returns flag of 0 if there has been no change, else 1
%returns the key that was used on exit in "button" 
%returns new left and right axis values

if nargin == 0
  do_patch = 1;
  h_axis_zoom = gca;
elseif nargin == 1
  do_patch = varargin{1};
  h_axis_zoom = gca;
elseif nargin > 1
  do_patch = varargin{1};
  h_axis_zoom = varargin{2};
end
  
%get new left and right values
[left, right, button, flag] = get_xrange(do_patch, h_axis_zoom);

%if there are changes, redisplay
if flag == 1;
  xlim=[left, right];
  set(gca,'xlim', xlim)
end
    
