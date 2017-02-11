 function [h_amp_plot, h_labels] = disp_song 
    %displays amplitude and label data in current plot window
    %displays data in the range set by axis limits saved with current axis handle
    % Now needs handles of amplitude and label axes since the are
    % separate now (BDW)
    
    %global variables for display
    
    global onsets offsets labels 
    global d_song time_d_song deriv_song
    global h_main_amp h_main_labels
    
    % if there are no onsets, then define one 20ms long syllable at 0.1 and
    % label it with '0' - Raghav 
    if (isempty(onsets))
        onsets = 100;
        offsets = 120;
        labels = '0';
    end
    
    %get value for scaling note display to song amplitude
    ymax=max(d_song);
    ymin = min(d_song);
    if (ymin > 0)
        ymin = 0;
    end
    
    %display     
    set(gcf, 'CurrentAxes', h_main_amp)
    cla
    hold on
    
    %plot song
    h_amp_plot = plot(time_d_song, d_song);
    set(h_amp_plot,'Tag','AmpPlot');
    
    %plot amplitude derivative - added by Raghav
    %h_amp_deriv_plot = plot(time_d_song, deriv_song/4 + 1.5, 'k');
    %set(h_amp_deriv_plot,'Tag','DerivPlot');
    
    %plot notes
    plot([onsets'; onsets'],[ymax*ones(size(onsets'));ymin*ones(size(onsets'))],'m');
    plot([offsets'; offsets'],[ymax*ones(size(offsets'));ymin*ones(size(offsets'))],'m');
    plot([onsets';offsets'],[ymax*ones(size(onsets'));ymax*ones(size(onsets'))],'m');

    hold off

    %plot labels
    
    set(gcf, 'CurrentAxes', h_main_labels)
    % Clear axes and then display new text 
    h_old_labels = findobj(h_main_labels,'Tag','SegLabels');
    if ( ~isempty(h_old_labels))
      delete(h_old_labels);
    end
    h_patch = findobj(h_main_labels,'Tag','LabelPatch');
    if ( ~isempty(h_patch))
      delete(h_patch);
    end

   
    %plot a white patch for labels to go on.
    left = min([time_d_song(1) onsets(1)])
    right  = max([time_d_song(end) offsets(end)])
    left = left - 10;
    right = right + 10;
    h_patch = patch('xdata',[left left right right],'ydata',[0.25 0.75 ...
                        0.75 0.25], 'FaceColor', 'w', 'LineStyle', ...
                    'none', 'Tag', 'LabelPatch');
    
    h_labels = text(onsets,ones(size(onsets))*0.5,labels','Clipping','on');
    set(h_labels,'Tag','SegLabels')
    
    set(gcf, 'CurrentAxes', h_main_amp)

    
