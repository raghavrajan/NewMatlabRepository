function uizoom(option, varargin)

%controls zooming and unzooming for songanal plots

if nargin == 1
  do_patch = 1;
  h_main = gcbf;
elseif nargin == 2
  do_patch = varargin{1};
  h_main = gcbf;
elseif nargin > 2
  do_patch = varargin{1};
  h_main = varargin{2};
end

SA = getappdata(h_main, 'SA_Data');
h_main_amp = SA.handles.mainAmp;
h_main_spect = SA.handles.mainSpect;
h_main_labels = SA.handles.mainLabels;

if strcmp(option,'x')
     %display x zoom message and zoom on amplitude
     
     subplot(h_main_amp);
     h_title = title('X-Zoom mode,select with mouse/ret, q to cancel','color',[0, 0, 1])
     [xmin, xmax, button, flag] = zoom_x(do_patch,h_main_amp);
     set(h_title,'String','');

     %if the amplitude plot was changed, also change the spectrogram plot
     %and vice-versa
     if flag == 1;
       h_curr = gca;
       ha = findobj(gcf,'Tag','SAAxis');
       for i=1:length(ha)
	 if ha(i) ~= h_curr
	   set(ha(i),'xlim', [xmin xmax])
	 end
       end
     end
elseif strcmp(option,'y')
       %display y zoom message and zoom on amplitude
       subplot(h_main_amp);
       title('Y-Zoom mode,select with mouse/ret, q to cancel','color',[0, 0, 1])
       [ymin, ymax, button, flag] = zoom_y(do_patch,h_main_amp);
       title('')     

elseif strcmp(option,'ux')
       %reset xmin & xmax, then redisplay 
       g_lim=get(h_main_amp,'userdata');
       set(h_main_amp,'xlim',[g_lim(1) g_lim(2)])
       set(h_main_labels,'xlim',[g_lim(1) g_lim(2)])
       g_lim=get(h_main_spect,'userdata');
       set(h_main_spect,'xlim',[g_lim(1) g_lim(2)])
       subplot(h_main_amp);

elseif strcmp(option,'uy')
       %reset ymin & ymax then redisplay
       subplot(h_main_amp);
       g_lim=get(gca,'userdata');
       set(gca,'ylim',[g_lim(3) g_lim(4)])    
end
