function gui_image(input1,input2,h,h_t,m_spec_deriv ,sound_bird,track,x,fs,  m_AM, m_FM ,m_Entropy , m_Amplitude , m_Freq , m_Pitch_gdns , m_Pitch ,  x_length , y_length  )

%version 1.0001
figure(h_t);
global cutoff  cutoff_value  axis_x1 axis_x2 axis_y1 axis_y2 add_record figure_table file_name bgcolor im_z fill_box save_location param
index_color_table=0;

bgcolor=[0.8 0.8 1];
cdata_index=round(adj_position([45,12]));
cdata_color=ones(cdata_index(2),cdata_index(1))*0.8;
cdata_color(:,:,1)=0.85;
cdata_color(:,:,2)=0.8;
cdata_color(:,:,3)=1;
ScreenSize=get(0,'ScreenSize');
        
v_matalb=ver('matlab');
if str2num(v_matalb.Version(1))~=7
    error('Matlab 7 is required to run ?sam?')
end
 
if isempty(x);

    load('Parameters');
    cutoff_value=param.cutoff_value;
    cutoff=param.cutoff;
    clear param
    
    track.duration_start=0;
    track.duration_end=0;
    track.duration_length=0;
    track.Pitch_mean=0;
    track.Pitch_min=0;
    track.Pitch_max=0;
    track.Pitch_var=0;
    track.Pitch_all=0;
    track.FM_mean=0;
    track.FM_min=0;
    track.FM_max=0;
    track.FM_var=0;
    track.FM_all=0;
    track.Entropy_mean=0;
    track.Entropy_min=0;
    track.Entropy_max=0;
    track.Entropy_var=0;
    track.Entropy_all=0;
    track.Pitch_gdns_mean=0;
    track.Pitch_gdns_min=0;
    track.Pitch_gdns_max=0;
    track.Pitch_gdns_var=0;
    track.Pitch_gdns_all=0;
    track.Amplitude_min=0;
    track.Amplitude_max=0;
    track.Amplitude_var=0;
    track.Amplitude_mean=0;
    track.Amplitude_all=0;
    track.M_Freq_min=0;
    track.M_Freq_max=0;
    track.M_Freq_var=0;
    track.M_Freq_mean=0;
    track.M_Freq_all=0;
    track.image=[];
    track.image_axes_x=[];
    track.image_axes_y=[];
    track.cutoff=[];
    fill_box=0;


        ScreenSize=get(0,'ScreenSize');
    %Color of line - ploting the features
        features_idn(1).color=[1 1 0];
        features_idn(1).position=120/768*ScreenSize(4);
        features_idn(1).str='Pitch';
        features_idn(2).color=[0 0.4 0];%[0.8 0 0.8];
        features_idn(2).position=100/768*ScreenSize(4);
        features_idn(2).str='FM';
        features_idn(3).color=[0 1 1];
        features_idn(3).position=80/768*ScreenSize(4);
        features_idn(3).str='Entropy';
        features_idn(4).color=[0 1 0];
        features_idn(4).position=60/768*ScreenSize(4);
        features_idn(4).str='Pitch gd';
        features_idn(5).color=[1 0.5 0.7];
        features_idn(5).position=40/768*ScreenSize(4);
        features_idn(5).str='Amplitude';
        features_idn(6).color=[1 0.5 0.1];
        features_idn(6).position=20/768*ScreenSize(4);
        features_idn(6).str='M Freq';



elseif nargin>7
    
    figure(h_t)
    h_line = findobj(gca,'Type','line');
    h_line_x=get(h_line,'XData');

    
    if im_z~=0
        delete(im_z)
        im_z=0;
    end
    uipanel('BackgroundColor',bgcolor,'Position',[0 .79 1 .21]);
    
    track.duration_start=x(1);%/fs*1000*44;
    track.duration_end=x(2);%/fs*1000*44;
    track.duration_length=(x(2)-x(1));%/fs*1000*44;
    
    load('Parameters');
    if param.pitch_averaging==0
        track.Pitch_mean=mean(m_Pitch(x(1):x(2)));%*44100/param.pad;
    else
        normelized_Pitch_gdns=m_Pitch(x(1):x(2)).*m_Pitch_gdns(x(1):x(2))/sum(m_Pitch_gdns(x(1):x(2)));%*44100/param.pad;
        track.Pitch_mean=sum(normelized_Pitch_gdns);
        
    end
    
        track.Pitch_var=var(m_Pitch(x(1):x(2)));%*44100/param.pad;    
    
    track.Pitch_min=min(m_Pitch(x(1):x(2)));%*44100/param.pad;
    track.Pitch_max=max(m_Pitch(x(1):x(2)));%*44100/param.pad;
    track.Pitch_all=(m_Pitch(x(1):x(2)));%*44100/param.pad;
    track.FM_mean=mean(m_FM(x(1):x(2)));
    track.FM_min=min(m_FM(x(1):x(2)));
    track.FM_max=max(m_FM(x(1):x(2)));
    track.FM_var=var(m_FM(x(1):x(2)));
    track.FM_all=(m_FM(x(1):x(2)));
    track.Entropy_mean=mean(m_Entropy(x(1):x(2)));
    track.Entropy_min=min(m_Entropy(x(1):x(2)));
    track.Entropy_max=max(m_Entropy(x(1):x(2)));
    track.Entropy_var=var(m_Entropy(x(1):x(2)));
    track.Entropy_all=(m_Entropy(x(1):x(2)));
    
    track.Pitch_gdns_mean=mean(m_Pitch_gdns(x(1):x(2)));
    track.Pitch_gdns_var=var(m_Pitch_gdns(x(1):x(2)));
    track.Pitch_gdns_min=min(m_Pitch_gdns(x(1):x(2)));
    track.Pitch_gdns_max=max(m_Pitch_gdns(x(1):x(2)));
    track.Pitch_gdns_all=(m_Pitch_gdns(x(1):x(2)));
    track.Amplitude_mean=mean(m_Amplitude(x(1):x(2)));
    track.Amplitude_min=min(m_Amplitude(x(1):x(2)));
    track.Amplitude_max=max(m_Amplitude(x(1):x(2)));
    track.Amplitude_var=var(m_Amplitude(x(1):x(2)));
    track.Amplitude_all=(m_Amplitude(x(1):x(2)));
    track.M_Freq_mean=mean(m_Freq(x(1):x(2)));
    track.M_Freq_min=min(m_Freq(x(1):x(2)));
    track.M_Freq_max=max(m_Freq(x(1):x(2)));
    track.M_Freq_var=var(m_Freq(x(1):x(2)));
    track.M_Freq_all=(m_Freq(x(1):x(2)));
    track.image=m_spec_deriv(:,x(1):x(2));
    track.image_axes_x=((x(1):x(2))*44.1);
    track.image_axes_y=(fs*(1:param.y_length)/param.pad);
    track.cutoff=cutoff_value;
    fill_box=1;

        ScreenSize=get(0,'ScreenSize');
    %Color of line - ploting the features
        features_idn(1).color=[1 1 0];
        features_idn(1).position=120/768*ScreenSize(4);
        features_idn(1).str='Pitch';
        features_idn(2).color=[0 0.4 0];%[0.8 0 0.8];
        features_idn(2).position=100/768*ScreenSize(4);
        features_idn(2).str='FM';
        features_idn(3).color=[0 1 1];
        features_idn(3).position=80/768*ScreenSize(4);
        features_idn(3).str='Entropy';
        features_idn(4).color=[0 1 0];
        features_idn(4).position=60/768*ScreenSize(4);
        features_idn(4).str='Pitch gd';
        features_idn(5).color=[1 0.5 0.7];
        features_idn(5).position=40/768*ScreenSize(4);
        features_idn(5).str='Amplitude';
        features_idn(6).color=[1 0.5 0.1];
        features_idn(6).position=20/768*ScreenSize(4);
        features_idn(6).str='M Freq';

index_color_table=[];
          figure(h);
            h_line = findobj(gca,'Type','line');
      if ~isempty(h_line)
         for index=1:6
            for k=1:length(h_line)
                cdata_line=get(h_line(k),'Color');
                if min(cdata_line==features_idn(index).color)
                    index_color_table=[index_color_table index];
                end
            end
         end
      end      


figure(h_t)

end
if 0 & ~isempty(h) 
    % ===The clear button===============================
    uicontrol('Style','pushbutton',...
        'Position',adj_position([222 1 45 15]),...
        'String','Clear',...
        'Interruptible','off',...
        'BusyAction','cancel',...
        'BackgroundColor',bgcolor,...
        'Callback',{@Clear,h,h_t,m_spec_deriv,[], x , fs , m_AM, m_FM ,m_Entropy , m_Amplitude ,m_Freq,  m_Pitch_gdns , m_Pitch   , x_length,y_length});
end
    % ===The play button===============================
    
    ration_sd=length(sound_bird)/max([size(m_spec_deriv) , 1]);

    
    play_song_image=imread('play1.jpg');
    play_song_image = imresiz(play_song_image,[size(play_song_image,1)*ScreenSize(4)/768 size(play_song_image,2)*ScreenSize(3)/1024]);
    if ~isempty(h)
    if ~isempty(x)
    uicontrol('Style','pushbutton',...
        'Position',adj_position([232 3 20 15]),...
        'Interruptible','off',...
        'CData',play_song_image,...
        'BusyAction','cancel',...
        'BackgroundColor',bgcolor,...
        'Callback',{@play_sound_bird,x,sound_bird,fs,ration_sd,h,h_t});
%        'TooltipString','Play the song of chosen segment',...
    else
    uicontrol('Style','pushbutton',...
        'Position',adj_position([232 3 20 15]),...
        'Interruptible','off',...
        'CData',play_song_image,...
        'BusyAction','cancel',...
        'BackgroundColor',bgcolor,...
        'Callback',{@play_sound_bird,x,sound_bird,fs,ration_sd,h,h_t});
%        'TooltipString','Play song',...
    end
    end
    
%        'String','Play',...
if ~isempty(h)

    % ===The Pan button===============================
    hand=imread('hand.jpg');
    hand = imresiz(hand,[size(hand,1)*ScreenSize(4)/768 size(hand,2)*ScreenSize(3)/1024]);

    uicontrol('Style','pushbutton',...
        'Position',adj_position([ 3 163 15 20]),...   %18 140 45 15]),...
        'Interruptible','off',...
        'CData',hand,...
        'TooltipString','Pan',...
        'BusyAction','cancel',...
        'BackgroundColor',bgcolor,...
        'Callback',{@Pan_x,h,h_t,bgcolor,m_spec_deriv,sound_bird, x_length , y_length , fs , m_AM, m_FM ,m_Entropy , m_Amplitude ,m_Freq,  m_Pitch_gdns , m_Pitch });
        %'String','Pan',...

end
if ~isempty(h)
    % ===The rhythm button===============================
    rhythm_sd=imread('rhythm.jpg');
    rhythm_sd = imresiz(rhythm_sd,[size(rhythm_sd,1)*ScreenSize(4)/768 size(rhythm_sd,2)*ScreenSize(3)/1024]);

    uicontrol('Style','pushbutton',...
        'Position',adj_position([30 141 20 16]),...
        'Interruptible','off',...
        'CData',rhythm_sd,...
           'TooltipString','User manual',...
        'BusyAction','cancel',...
        'BackgroundColor',bgcolor,...
        'Callback',{@rhythm_spectral_derivatives,fs,  m_AM, m_FM ,m_Entropy , m_Amplitude , m_Freq , m_Pitch_gdns , m_Pitch });
        %'String','Pan',...
end
    % ===Store record button===============================
    open_song_image=imread('open.jpg');
    open_song_image = imresiz(open_song_image,[size(open_song_image,1)*ScreenSize(4)/768 size(open_song_image,2)*ScreenSize(3)/1024]);
     if ~isempty(h)
    save_image=imread('save.jpg');
    save_image = imresiz(save_image,[size(save_image,1)*ScreenSize(4)/768 size(save_image,2)*ScreenSize(3)/1024]);

         uicontrol('Style','pushbutton',...
        'Position',adj_position([132 193 16 15]),...
        'Interruptible','off',...
        'BusyAction','cancel',...
        'TooltipString','Store information of the chosen segment',... 
        'CData',save_image,...
         'BackgroundColor',bgcolor,...
        'Callback',{@Add_record,track,h});
    %'String','Store',...
        
     end
if 0
     % ===Axis button===============================
     if ~isempty(h)
    uicontrol('Style','pushbutton',...
        'String','Axis',...
        'Position',adj_position([169 193 44 15]),...
        'Interruptible','off',...
        'BusyAction','cancel',...
        'BackgroundColor',bgcolor,...
        'Callback',{@plot_axis,h,h_t,m_spec_deriv});
     end
end
     % ===Attach button===============================
     if ~isempty(h)
    attach_song_image=imread('attach.jpg');
    attach_song_image = imresiz(attach_song_image,[size(attach_song_image,1)*ScreenSize(4)/768 size(attach_song_image,2)*ScreenSize(3)/1024]);
    uicontrol('Style','pushbutton',...
        'Position',adj_position([177 193 30 12]),...
        'Interruptible','off',...
        'CData',attach_song_image,...
        'TooltipString','Restore original windows position of current figure',... 
        'BusyAction','cancel',...
        'BackgroundColor',bgcolor,...
        'Callback',{@Attach,h,h_t});
%            'String','Align',...
     end
    % ===Parameters button===============================
    P_song_image=imread('P.jpg');
    P_song_image = imresiz(P_song_image,[size(P_song_image,1)*ScreenSize(4)/768 size(P_song_image,2)*ScreenSize(3)/1024]);
     
    uicontrol('Style','pushbutton',...
        'Position',adj_position([237 193 15 17]),...
        'CData',P_song_image,...
        'Interruptible','off',...
        'TooltipString','Change parameters',... 
        'BusyAction','cancel',...
        'BackgroundColor',bgcolor,...
        'Callback',{@Change_Parameters,44100});
 %      'String','P',...
 
    % ===Open button===============================
    open_song_image=imread('open.jpg');
    open_song_image = imresiz(open_song_image,[size(open_song_image,1)*ScreenSize(4)/768 size(open_song_image,2)*ScreenSize(3)/1024]);
     
    uicontrol('Style','pushbutton',...
        'Position',adj_position([32 193 15 15]),...
        'Interruptible','off',...
        'BusyAction','cancel',...
        'TooltipString','Open file in a new window',... 
           'CData',open_song_image,...
        'BackgroundColor',bgcolor,...
        'Callback',{@Open_song,h,h_t});
    %'String','Open',...
        
     if ~isempty(h)

         position_index=getappdata(h_t,'position_index');

         
    re_open_song_image=imread('reopen4.jpg');
    re_open_song_image = imresiz(re_open_song_image,[size(re_open_song_image,1)*ScreenSize(4)/768 size(re_open_song_image,2)*ScreenSize(3)/1024]);

        uicontrol('Style','pushbutton',...
        'Position',adj_position([82 192 18 15]),...
        'Interruptible','off',...
        'TooltipString','Open file in the same window',... 
        'BusyAction','cancel',...
           'CData',re_open_song_image,...
        'BackgroundColor',bgcolor,...
        'Callback',{@Open_song,h,h_t,position_index});
     end

    if 0
    % ===The quit button===============================
     
    uicontrol('Style','pushbutton',...
        'Position',adj_position([18 140 45 15]),...
        'String','Quit',...
        'Interruptible','off',...
        'BusyAction','cancel',...
        'BackgroundColor',bgcolor,...
        'Callback',@quit_program);
   end
    if ~isempty(h)
   % ===The contrast slider===============================
     
    uicontrol('Style','slider',...
        'Position',adj_position([1 20 11 135]),... %755= 0
        'Min',1,'Max',10,'value',cutoff_value,...
        'SliderStep',[0.01 0.1],...
        'String','Contrast',...
        'Interruptible','off',...
        'BusyAction','cancel',...
        'BackgroundColor',[0.4 0.4 0.8],...
        'Callback',{@contrast_adj,h,h_t});
    end

    % ===image==========
axes('Units','pixel','position', adj_position([17 163  257 19]))
h_i=imagesc(m_spec_deriv);   axis off;
   colormap(bone); set(gca,'YDir','normal'); caxis([-cutoff cutoff])
if ~isempty(x)
hold on
plot(h_line_x,[0 500],'r');
end
if fill_box
    uicontrol('Style','text',...
        'Position',adj_position([30 3 200 15]),...
        'String','                                     ',...
    'BackgroundColor',bgcolor);

    % ===text==========
    uicontrol(gcf,'Style','text',...
        'String', 'Mean',...
        'Position',adj_position([65 138 50 18]),...
        'BackgroundColor',bgcolor);
    uicontrol(gcf,'Style','text',...
        'String', 'Min',...
        'Position',adj_position([115 138 50 18]),...
        'BackgroundColor',bgcolor);
    uicontrol(gcf,'Style','text',...
        'String', 'Max',...
        'Position',adj_position([165 138 50 18]),...
        'BackgroundColor',bgcolor);
    uicontrol(gcf,'Style','text',...
        'String', 'Var',...
        'Position',adj_position([215 138 50 18]),...
        'BackgroundColor',bgcolor);

    uicontrol(gcf,'Style','text',...
        'Position',adj_position([265 138 50 18]),...
        'BackgroundColor',bgcolor);
end
if ~isempty(h)

        
    if 0 %isempty(h)
        uicontrol(gcf,'Style','text',...
            'String', 'Pitch',...
            'Position',adj_position([15 120 50 18]),...
            'BackgroundColor',bgcolor);
        uicontrol(gcf,'Style','text',...
            'String', 'FM',...
            'Position',adj_position([15 100 50 18]),...
            'BackgroundColor',bgcolor);
        uicontrol(gcf,'Style','text',...
            'String', 'Entropy',...
            'Position',adj_position([15 80 50 18]),...
            'BackgroundColor',bgcolor);
        uicontrol(gcf,'Style','text',...
            'String', 'Pitch gdns',...
            'Position',adj_position([15 60 50 18]),...
            'BackgroundColor',bgcolor);
        uicontrol(gcf,'Style','text',...
            'String', 'Amplitude',...
            'Position',adj_position([15 40 50 18]),...
            'BackgroundColor',bgcolor);
        uicontrol(gcf,'Style','text',...
            'String', 'M Freq',...
            'Position',adj_position([15 20 50 18]),...
            'BackgroundColor',bgcolor);
    end        
        if max(index_color_table==1)
         index=1;
         tmp_cdata_color=cdata_color;
         tmp_cdata_color(:,:,1)=features_idn(index).color(1);
         tmp_cdata_color(:,:,2)=features_idn(index).color(2);
         tmp_cdata_color(:,:,3)=features_idn(index).color(3);
        uicontrol(gcf,'Style','pushbutton',...
            'String', 'Pitch',...
            'Interruptible','off',...
            'CData',tmp_cdata_color,...
            'BusyAction','cancel',...
            'Callback',{@plot_feature,x_length,y_length,m_Pitch,features_idn,1,h,h_t,cdata_color },...
            'Position',adj_position([15 120 50 18]),...
            'BackgroundColor',bgcolor);
        else
                    uicontrol(gcf,'Style','pushbutton',...
            'String', 'Pitch',...
            'Interruptible','off',...
            'CData',cdata_color,...
            'BusyAction','cancel',...
            'Callback',{@plot_feature,x_length,y_length,m_Pitch,features_idn,1,h,h_t,cdata_color },...
            'Position',adj_position([15 120 50 18]),...
            'BackgroundColor',bgcolor);
        end
        
        if max(index_color_table==2)
        index=2;
         tmp_cdata_color=cdata_color;
         tmp_cdata_color(:,:,1)=features_idn(index).color(1);
         tmp_cdata_color(:,:,2)=features_idn(index).color(2);
         tmp_cdata_color(:,:,3)=features_idn(index).color(3);
        uicontrol(gcf,'Style','pushbutton',...
            'String', 'FM',...
            'Interruptible','off',...
            'CData',tmp_cdata_color,...
            'BusyAction','cancel',...
            'Callback',{@plot_feature,x_length,y_length,m_FM,features_idn,2,h,h_t,cdata_color},...
            'Position',adj_position([15 100 50 18]),...
            'BackgroundColor',bgcolor);
        else
        uicontrol(gcf,'Style','pushbutton',...
            'String', 'FM',...
            'Interruptible','off',...
            'CData',cdata_color,...
            'BusyAction','cancel',...
            'Callback',{@plot_feature,x_length,y_length,m_FM,features_idn,2,h,h_t,cdata_color},...
            'Position',adj_position([15 100 50 18]),...
            'BackgroundColor',bgcolor);
        end
        if max(index_color_table==3)
        index=3;
         tmp_cdata_color=cdata_color;
         tmp_cdata_color(:,:,1)=features_idn(index).color(1);
         tmp_cdata_color(:,:,2)=features_idn(index).color(2);
         tmp_cdata_color(:,:,3)=features_idn(index).color(3);
        uicontrol(gcf,'Style','pushbutton',...
            'String', 'Entropy',...
            'Interruptible','off',...
            'CData',tmp_cdata_color,...
            'BusyAction','cancel',...
            'Callback',{@plot_feature,x_length,y_length,m_Entropy,features_idn,3,h,h_t,cdata_color},...
            'Position',adj_position([15 80 50 18]),...
            'BackgroundColor',bgcolor);
        else
            uicontrol(gcf,'Style','pushbutton',...
            'String', 'Entropy',...
            'Interruptible','off',...
            'CData',cdata_color,...
            'BusyAction','cancel',...
            'Callback',{@plot_feature,x_length,y_length,m_Entropy,features_idn,3,h,h_t,cdata_color},...
            'Position',adj_position([15 80 50 18]),...
            'BackgroundColor',bgcolor);
        end
        if max(index_color_table==4)
        index=4;
         tmp_cdata_color=cdata_color;
         tmp_cdata_color(:,:,1)=features_idn(index).color(1);
         tmp_cdata_color(:,:,2)=features_idn(index).color(2);
         tmp_cdata_color(:,:,3)=features_idn(index).color(3);
        uicontrol(gcf,'Style','pushbutton',...
            'String', 'Pitch gd',...
            'Interruptible','off',...
            'BusyAction','cancel',...
            'CData',tmp_cdata_color,...
            'Callback',{@plot_feature,x_length,y_length,m_Pitch_gdns,features_idn,4,h,h_t,cdata_color},...
            'Position',adj_position([15 60 50 18]),...
            'BackgroundColor',bgcolor);
        else
            uicontrol(gcf,'Style','pushbutton',...
            'String', 'Pitch gd',...
            'Interruptible','off',...
            'BusyAction','cancel',...
            'CData',cdata_color,...
            'Callback',{@plot_feature,x_length,y_length,m_Pitch_gdns,features_idn,4,h,h_t,cdata_color},...
            'Position',adj_position([15 60 50 18]),...
            'BackgroundColor',bgcolor);
        end
        if max(index_color_table==5)
        index=5;
         tmp_cdata_color=cdata_color;
         tmp_cdata_color(:,:,1)=features_idn(index).color(1);
         tmp_cdata_color(:,:,2)=features_idn(index).color(2);
         tmp_cdata_color(:,:,3)=features_idn(index).color(3);
        uicontrol(gcf,'Style','pushbutton',...
            'String', 'Amplitude',...
            'Interruptible','off',...
            'BusyAction','cancel',...
            'Callback',{@plot_feature,x_length,y_length,m_Amplitude,features_idn,5,h,h_t,cdata_color},...
            'Position',adj_position([15 40 50 18]),...
            'CData',tmp_cdata_color,...
            'BackgroundColor',bgcolor);
        else
            uicontrol(gcf,'Style','pushbutton',...
            'String', 'Amplitude',...
            'Interruptible','off',...
            'BusyAction','cancel',...
            'Callback',{@plot_feature,x_length,y_length,m_Amplitude,features_idn,5,h,h_t,cdata_color},...
            'Position',adj_position([15 40 50 18]),...
            'CData',cdata_color,...
            'BackgroundColor',bgcolor);
        end
        if max(index_color_table==6)
        index=6;
         tmp_cdata_color=cdata_color;
         tmp_cdata_color(:,:,1)=features_idn(index).color(1);
         tmp_cdata_color(:,:,2)=features_idn(index).color(2);
         tmp_cdata_color(:,:,3)=features_idn(index).color(3);
        uicontrol(gcf,'Style','pushbutton',...
            'String', 'M Freq',...
            'Interruptible','off',...
            'BusyAction','cancel',...
            'Callback',{@plot_feature,x_length,y_length,m_Freq,features_idn,6,h,h_t,cdata_color},...
            'Position',adj_position([15 20 50 18]),...
            'CData',tmp_cdata_color,...
            'BackgroundColor',bgcolor);
        else
            uicontrol(gcf,'Style','pushbutton',...
            'String', 'M Freq',...
            'Interruptible','off',...
            'BusyAction','cancel',...
            'Callback',{@plot_feature,x_length,y_length,m_Freq,features_idn,6,h,h_t,cdata_color},...
            'Position',adj_position([15 20 50 18]),...
            'CData',cdata_color,...
            'BackgroundColor',bgcolor);
        end
    if ~isempty(x)
    
    uicontrol(gcf,'Style','text',...
        'String', 'Duration',...
        'Position',adj_position([15 1 50 18]),...
        'BackgroundColor',bgcolor);
    end
end
if fill_box

%=================features
    % ===Pitch============================================
    uicontrol(gcf,'Style','edit',...
        'String', sprintf( '%.0f ', track.Pitch_mean ),...
        'Position',adj_position([70 123 47 14]),...
        'BackgroundColor',[1 1 1],...         
        'CallBack', {@gui_image,h,h_t,m_spec_deriv,sound_bird,[], x , fs , m_AM, m_FM ,m_Entropy , m_Amplitude ,m_Freq,  m_Pitch_gdns , m_Pitch ,x_length,y_length }); 
    

    % ===FM============================================
    uicontrol(gcf,'Style','edit',...
        'String', sprintf( '%.1f ', track.FM_mean ),...
        'Position',adj_position([70 103 47 14]),...
        'BackgroundColor',[1 1 1],...         
        'CallBack', {@gui_image,h,h_t,m_spec_deriv,sound_bird,[], x , fs , m_AM, m_FM ,m_Entropy , m_Amplitude ,m_Freq,  m_Pitch_gdns , m_Pitch  ,x_length,y_length }); 

    % ===Entropy============================================
    uicontrol(gcf,'Style','edit',...
        'String', sprintf( '%.1f ', track.Entropy_mean ),...
        'Position',adj_position([70 83 47 14]),...
        'BackgroundColor',[1 1 1],...         
        'CallBack', {@gui_image,h,h_t,m_spec_deriv,sound_bird,[], x , fs , m_AM, m_FM ,m_Entropy , m_Amplitude ,m_Freq,  m_Pitch_gdns , m_Pitch  ,x_length,y_length }); 

    % ===Pitch_gdns============================================
    
    
    uicontrol(gcf,'Style','edit',...
        'String', sprintf( '%.0f ', track.Pitch_gdns_mean ),...
        'Position',adj_position([70 63 47 14]),...
        'BackgroundColor',[1 1 1],...         
        'CallBack', {@gui_image,h,h_t,m_spec_deriv,sound_bird,[], x , fs , m_AM, m_FM ,m_Entropy , m_Amplitude ,m_Freq,  m_Pitch_gdns , m_Pitch  ,x_length,y_length }); 


    % ===Amplitude============================================
    uicontrol(gcf,'Style','edit',...
        'String', sprintf( '%.1f ', track.Amplitude_mean ),...
        'Position',adj_position([70 43 47 14]),...
        'BackgroundColor',[1 1 1],...
        'CallBack', {@gui_image,h,h_t,m_spec_deriv,sound_bird,[], x , fs , m_AM, m_FM ,m_Entropy , m_Amplitude ,m_Freq,  m_Pitch_gdns , m_Pitch  ,x_length,y_length }); 


    % ===Freq============================================
    uicontrol(gcf,'Style','edit',...  
        'String', sprintf( '%.0f ', track.M_Freq_mean ),...
        'Position',adj_position([70 23 47 14]),...
        'BackgroundColor',[1 1 1],...
        'CallBack', {@gui_image,h,h_t,m_spec_deriv,sound_bird,[], x , fs , m_AM, m_FM ,m_Entropy , m_Amplitude ,m_Freq,  m_Pitch_gdns , m_Pitch  ,x_length,y_length }); 
%==============================


    % ===Pitch============================================
    uicontrol(gcf,'Style','edit',...
        'String', sprintf( '%.0f ', track.Pitch_min ),...
        'Position',adj_position([120 123 47 14]),...
        'BackgroundColor',[1 1 1],...
        'CallBack', {@gui_image,h,h_t,m_spec_deriv,sound_bird,[], x , fs , m_AM, m_FM ,m_Entropy , m_Amplitude ,m_Freq,  m_Pitch_gdns , m_Pitch  ,x_length,y_length }); 
    

    % ===FM============================================
    uicontrol(gcf,'Style','edit',...
        'String', sprintf( '%.1f ', track.FM_min ),...
        'Position',adj_position([120 103 47 14]),...
        'BackgroundColor',[1 1 1],...
        'CallBack', {@gui_image,h,h_t,m_spec_deriv,sound_bird,[], x , fs , m_AM, m_FM ,m_Entropy , m_Amplitude ,m_Freq,  m_Pitch_gdns , m_Pitch  ,x_length,y_length }); 

    % ===Entropy============================================
    uicontrol(gcf,'Style','edit',...
        'String', sprintf( '%.1f ', track.Entropy_min ),...
        'Position',adj_position([120 83 47 14]),...
        'BackgroundColor',[1 1 1],...
        'CallBack', {@gui_image,h,h_t,m_spec_deriv,sound_bird,[], x , fs , m_AM, m_FM ,m_Entropy , m_Amplitude ,m_Freq,  m_Pitch_gdns , m_Pitch  ,x_length,y_length }); 

    % ===Pitch_gdns============================================
    uicontrol(gcf,'Style','edit',...
        'String', sprintf( '%.0f ', track.Pitch_gdns_min ),...
        'Position',adj_position([120 63 47 14]),...
        'BackgroundColor',[1 1 1],...
        'CallBack', {@gui_image,h,h_t,m_spec_deriv,sound_bird,[], x , fs , m_AM, m_FM ,m_Entropy , m_Amplitude ,m_Freq,  m_Pitch_gdns , m_Pitch  ,x_length,y_length }); 


    % ===Amplitude============================================
    uicontrol(gcf,'Style','edit',...
        'String', sprintf( '%.1f ', track.Amplitude_min ),...
        'Position',adj_position([120 43 47 14]),...
        'BackgroundColor',[1 1 1],...
        'CallBack', {@gui_image,h,h_t,m_spec_deriv,sound_bird,[], x , fs , m_AM, m_FM ,m_Entropy , m_Amplitude ,m_Freq,  m_Pitch_gdns , m_Pitch  ,x_length,y_length }); 


    % ===Freq============================================
    uicontrol(gcf,'Style','edit',...
        'String', sprintf( ' %.0f', track.M_Freq_min ),...
        'Position',adj_position([120 23 47 14]),...
        'BackgroundColor',[1 1 1],...
        'CallBack', {@gui_image,h,h_t,m_spec_deriv,sound_bird,[], x , fs , m_AM, m_FM ,m_Entropy , m_Amplitude ,m_Freq,  m_Pitch_gdns , m_Pitch  ,x_length,y_length }); 
%======================================


    % ===Pitch============================================
    uicontrol(gcf,'Style','edit',...
        'String', sprintf( '%.0f ', track.Pitch_max ),...
        'Position',adj_position([170 123 47 14]),...
        'BackgroundColor',[1 1 1],...
        'CallBack', {@gui_image,h,h_t,m_spec_deriv,sound_bird,[], x , fs , m_AM, m_FM ,m_Entropy , m_Amplitude ,m_Freq,  m_Pitch_gdns , m_Pitch  ,x_length,y_length }); 
    

    % ===FM============================================
    uicontrol(gcf,'Style','edit',...
        'String', sprintf( '%.1f ', track.FM_max ),...
        'Position',adj_position([170 103 47 14]),...
        'BackgroundColor',[1 1 1],...
        'CallBack', {@gui_image,h,h_t,m_spec_deriv,sound_bird,[], x , fs , m_AM, m_FM ,m_Entropy , m_Amplitude ,m_Freq,  m_Pitch_gdns , m_Pitch  ,x_length,y_length }); 

    % ===Entropy============================================
    uicontrol(gcf,'Style','edit',...
        'String', sprintf( '%.1f ', track.Entropy_max ),...
        'Position',adj_position([170 83 47 14]),...
        'BackgroundColor',[1 1 1],...
        'CallBack', {@gui_image,h,h_t,m_spec_deriv,sound_bird,[], x , fs , m_AM, m_FM ,m_Entropy , m_Amplitude ,m_Freq,  m_Pitch_gdns , m_Pitch  ,x_length,y_length }); 

    % ===Pitch_gdns============================================
    uicontrol(gcf,'Style','edit',...
        'String', sprintf( '%.0f ', track.Pitch_gdns_max ),...
        'Position',adj_position([170 63 47 14]),...
        'BackgroundColor',[1 1 1],...
        'CallBack', {@gui_image,h,h_t,m_spec_deriv,sound_bird,[], x , fs , m_AM, m_FM ,m_Entropy , m_Amplitude ,m_Freq,  m_Pitch_gdns , m_Pitch  ,x_length,y_length }); 


    % ===Amplitude============================================
    uicontrol(gcf,'Style','edit',...
        'String', sprintf( '%.1f ', track.Amplitude_max ),...
        'Position',adj_position([170 43 47 14]),...
        'BackgroundColor',[1 1 1],...
        'CallBack', {@gui_image,h,h_t,m_spec_deriv,sound_bird,[], x , fs , m_AM, m_FM ,m_Entropy , m_Amplitude ,m_Freq,  m_Pitch_gdns , m_Pitch  ,x_length,y_length }); 


    % ===Freq============================================
    uicontrol(gcf,'Style','edit',...
        'String', sprintf( '%.0f ', track.M_Freq_max ),...
        'Position',adj_position([170 23 47 14]),...
        'BackgroundColor',[1 1 1],...
        'CallBack', {@gui_image,h,h_t,m_spec_deriv,sound_bird,[], x , fs , m_AM, m_FM ,m_Entropy , m_Amplitude ,m_Freq,  m_Pitch_gdns , m_Pitch  ,x_length,y_length }); 
%==================================


    % ===Pitch============================================
    uicontrol(gcf,'Style','edit',...
        'String', sprintf( '%.0f ', track.Pitch_var ),...
        'Position',adj_position([220 123 47 14]),...
        'TooltipString',num2str(track.Pitch_var),...
        'BackgroundColor',[1 1 1],...
        'CallBack', {@gui_image,h,h_t,m_spec_deriv,sound_bird,[], x , fs , m_AM, m_FM ,m_Entropy , m_Amplitude ,m_Freq,  m_Pitch_gdns , m_Pitch  ,x_length,y_length }); 
        

    % ===FM============================================
    uicontrol(gcf,'Style','edit',...
        'String', sprintf( '%.1f ', track.FM_var ),...
        'Position',adj_position([220 103 47 14]),...
        'TooltipString',num2str(track.FM_var),...
        'BackgroundColor',[1 1 1],...
        'CallBack', {@gui_image,h,h_t,m_spec_deriv,sound_bird,[], x , fs , m_AM, m_FM ,m_Entropy , m_Amplitude ,m_Freq,  m_Pitch_gdns , m_Pitch  ,x_length,y_length }); 
        
    % ===Entropy============================================
    uicontrol(gcf,'Style','edit',...
        'String', sprintf( '%.1f ', track.Entropy_var ),...
        'Position',adj_position([220 83 47 14]),...
        'TooltipString',num2str(track.Entropy_var),...
        'BackgroundColor',[1 1 1],...
        'CallBack', {@gui_image,h,h_t,m_spec_deriv,sound_bird,[], x , fs , m_AM, m_FM ,m_Entropy , m_Amplitude ,m_Freq,  m_Pitch_gdns , m_Pitch  ,x_length,y_length }); 
        
    % ===Pitch_gdns============================================
    uicontrol(gcf,'Style','edit',...
        'String', sprintf( '%.0f', track.Pitch_gdns_var ),...
        'Position',adj_position([220 63 47 14]),...
        'TooltipString',num2str(track.Pitch_gdns_var),...
        'BackgroundColor',[1 1 1],...
        'CallBack', {@gui_image,h,h_t,m_spec_deriv,sound_bird,[], x , fs , m_AM, m_FM ,m_Entropy , m_Amplitude ,m_Freq,  m_Pitch_gdns , m_Pitch  ,x_length,y_length }); 
        

    % ===Amplitude============================================
    uicontrol(gcf,'Style','edit',...
        'String', sprintf( '%.1f ', track.Amplitude_var ),...
        'Position',adj_position([220 43 47 14]),...
        'BackgroundColor',[1 1 1],...  
        'TooltipString',num2str(track.Amplitude_var),...
        'CallBack', {@gui_image,h,h_t,m_spec_deriv,sound_bird,[], x , fs , m_AM, m_FM ,m_Entropy , m_Amplitude ,m_Freq,  m_Pitch_gdns , m_Pitch  ,x_length,y_length }); 
        

    % ===Freq============================================
    uicontrol(gcf,'Style','edit',...
        'String', sprintf( '%.0f ', track.M_Freq_var ),...
        'Position',adj_position([220 23 47 14]),...
        'TooltipString',num2str(track.M_Freq_var),...
        'BackgroundColor',[1 1 1],... 
        'CallBack', {@gui_image,h,h_t,m_spec_deriv,sound_bird,[], x , fs , m_AM, m_FM ,m_Entropy , m_Amplitude ,m_Freq,  m_Pitch_gdns , m_Pitch  ,x_length,y_length }); 
        
    
        
        
        
        
    % ===Duration============================================
    uicontrol(gcf,'Style','edit',...
        'String', sprintf( '%.1f ', track.duration_length ),...
        'Position',adj_position([70 3 47 14]),...
        'BackgroundColor',[1 1 1],...
        'CallBack', {@gui_image,h,h_t,m_spec_deriv ,track });

    uicontrol(gcf,'Style','edit',...
        'String', sprintf( '%.1f ', track.duration_start ),...
        'Position',adj_position([120 3 47 14]),...
        'BackgroundColor',[1 1 1],...
        'CallBack', {@gui_image,h,h_t,m_spec_deriv ,track });

    uicontrol(gcf,'Style','edit',...
        'String', sprintf( '%.1f ', track.duration_end ),...
        'Position',adj_position([170 3 47 14]),...
        'BackgroundColor',[1 1 1],...
        'CallBack', {@gui_image,h,h_t,m_spec_deriv ,track });

elseif isempty(h)
    zebra_finch_carton=imread('zebra_finch_carton_plastic.jpg');
    zebra_finch_carton = imresiz(zebra_finch_carton,[size(zebra_finch_carton,1)*ScreenSize(4)/768 size(zebra_finch_carton,2)*ScreenSize(3)/1024]);
        im_z=uicontrol('Style','pushbutton',...
            'Position',adj_position([75 20 170 130]),...
              'Interruptible','off',...
            'BusyAction','cancel',...
            'CData',zebra_finch_carton);
        
    zebra_finch_carton=imread('sam.jpg');
    zebra_finch_carton = imresiz(zebra_finch_carton,[size(zebra_finch_carton,1)*ScreenSize(4)/768 size(zebra_finch_carton,2)*ScreenSize(3)/1024]);
        uicontrol('Style','pushbutton',...
            'Position',adj_position([35 3 176 14]),...
              'Interruptible','off',...
            'BusyAction','cancel',...
            'CData',zebra_finch_carton);
            if 0
                uicontrol('Style','text',...
                    'Position',adj_position([30 3 180 14]),...
                'FontAngle','normal',...
                'FontWeight','bold',...
                'String','Sound Anlaysis based on Matlab',...
                'BackgroundColor',bgcolor);
            end
        
else
    
    zebra_finch_carton=imread('credit.jpg');
    zebra_finch_carton = imresiz(zebra_finch_carton,[size(zebra_finch_carton,1)*ScreenSize(4)/768 size(zebra_finch_carton,2)*ScreenSize(3)/1024]);
        im_z=uicontrol('Style','pushbutton',...
            'Position',adj_position([65 20 200 134]),...%([85 20 170 130]),...
              'Interruptible','off',...
            'BusyAction','cancel',...
            'CData',zebra_finch_carton);
        
    zebra_finch_carton=imread('sam.jpg');
    zebra_finch_carton = imresiz(zebra_finch_carton,[size(zebra_finch_carton,1)*ScreenSize(4)/768 size(zebra_finch_carton,2)*ScreenSize(3)/1024]);
        uicontrol('Style','pushbutton',...
            'Position',adj_position([15 3 210 14]),...
              'Interruptible','off',...
            'BusyAction','cancel',...
            'CData',zebra_finch_carton);
            if 0
                uicontrol('Style','text',...
                    'Position',adj_position([30 3 100 15]),...
                'FontAngle','normal',...
                'FontWeight','bold',...
                'String','Sound Anlaysis based on Matlab',...
                'BackgroundColor',bgcolor);
            end
    
end
    
        
    %===== Box at the top figure ==============
    figure(h_t)
    uipanel('BackgroundColor',bgcolor,'Position',[0 .79 1 .21]);
        
   
   if nargin==0
        
        add_record.duration_start=[];
        add_record.duration_end=[];
        add_record.duration_length=[];
        add_record.Pitch_mean=[];
        add_record.Pitch_min=[];
        add_record.Pitch_max=[];
        add_record.Pitch_var=[];
        add_record.FM_mean=[];
        add_record.FM_min=[];
        add_record.FM_max=[];
        add_record.FM_var=[];
        add_record.Entropy_mean=[];
        add_record.Entropy_min=[];
        add_record.Entropy_max=[];
        add_record.Entropy_var=[];
        add_record.Pitch_gdns_mean=[];
        add_record.Pitch_gdns_min=[];
        add_record.Pitch_gdns_max=[];
        add_record.Pitch_gdns_var=[];
        add_record.Amplitude_min=[];
        add_record.Amplitude_max=[];
        add_record.Amplitude_var=[];
        add_record.Amplitude_mean=[];
        add_record.M_Freq_min=[];
        add_record.M_Freq_max=[];
        add_record.M_Freq_var=[];
        add_record.M_Freq_mean=[];
        add_record.image=[];
        add_record.image_axes_x=[];
        add_record.image_axes_y=[];
        add_record.comment=[];
   end
    
 

%=========Contrast slider callback===========


function []=contrast_adj(hObject,eventdata,h,h_t)

global cutoff  cutoff_value

load('Parameters');

figure(h_t);
h_tt=gcbo;

change_contrast=get(h_tt,'value');
cutoff_value=change_contrast;
if change_contrast<param.cutoff_value
    change_contrast=exp((param.cutoff_value-change_contrast));
    cutoff=50*(change_contrast/2+.5)*758400000;
else
%    change_contrast=1-change_contrast;
%    change_contrast=(log((change_contrast-param.cutoff_value+1)));
    change_contrast=1/exp((change_contrast-param.cutoff_value));
%    change_contrast=    change_contrast/max(change_contrast);
    cutoff=50*(change_contrast/2)*758400000;
end

figure(h);
caxis([-cutoff cutoff]);

return

%=========Quit callback===========


function quit_program(hObject,eventdata)
    global quit_ntn
quit_ntn=0;

return

%=========Clear callback===========


function []=Clear(hObject,eventdata,h,h_t,m_spec_deriv,track, x , fs , m_AM, m_FM ,m_Entropy , m_Amplitude ,m_Freq,  m_Pitch_gdns , m_Pitch  ,x_length,y_length)

global cutoff bgcolor
%cdata_index=round(adj_position([45,15]));
%cdata_color=ones(cdata_index(2),cdata_index(1))*0.8;
%cdata_color(:,:,1)=0.7;
%cdata_color(:,:,2)=0.8;
%cdata_color(:,:,3)=0.4;
        
        figure(h_t)
        
        h_line = findobj(gca,'Type','line');
        h_ui = findobj(gca,'Type','uicontrol');
        pan_on=0;
        for k=1:length(h_ui)
            if ~min(min(get(h_ui(k),'CData')==cdata_color))
                pan_on=1
            end
        end

      if ~isempty(h_line) 
        save_line=get(h_line);

        gui_image([],[],h,h_t,m_spec_deriv,sound_bird,[], x , fs , m_AM, m_FM ,m_Entropy , m_Amplitude ,m_Freq,  m_Pitch_gdns , m_Pitch  ,x_length,y_length);

        hold on
        plot(save_line.XData,save_line.YData,'Color',save_line.Color)
        if pan_on
             hand=imread('hand_g.jpg');
            hand = imresiz(hand,[size(hand,1)*ScreenSize(4)/768 size(hand,2)*ScreenSize(3)/1024]);

            uicontrol('Style','pushbutton',...
                'Position',adj_position([ 1 163 17 20]),...   %18 140 45 15]),...
                'Interruptible','off',...
                'CData',hand,...
                'BusyAction','cancel',...
                'BackgroundColor',bgcolor,...
                'Callback',{@Pan_x,h,h_t,bgcolor});
        end

            
       end
          figure(h);
            h_line = findobj(gca,'Type','line');
            for k=1:length(h_line)
                xdata=get(h_line(k),'XData');
                if length(xdata)>2
                    delete(h_line(k));
                end
            end
            
return

%=========Add record callback===========


function []=Add_record(hObject,eventdata,track,h)

global add_record save_location file_name store_path save_location figure_table

    if isempty(save_location)
        
        prompt={'Create a new structure? y/n '};
        name='Add record?';
        numlines=1;
        defaultanswer={'Yes'};

        answer=inputdlg(prompt,name,numlines,defaultanswer);
        options.Resize='on';
        options.WindowStyle='normal';
        options.Interpreter='tex';

        if ~isempty(answer)
            tmp_struct=cell2struct(answer(1),'tmp',1); %save
            New_struct=(tmp_struct.tmp(1));

            if New_struct=='Y' | New_struct=='y'
                curr_cd=cd;
                store_path = uigetdir(curr_cd, 'Pick a Directory to save in');
                save_location=[store_path '\Store_record.mat']; 
            else
                curr_cd=cd;
                
                [matname, mat_pathname, filterindex] = uigetfile( ...
                   {'*.mat','MAT-files (*.mat)'; ...
                    '*.*',  'All Files (*.*)'}, ...
                    'Pick saved structure file', ...
                    'MultiSelect', 'on');
                save_location=[mat_pathname matname]; 
            end
    
        end
    end  
    if     store_path~=0
        prompt={'Store location','Add Comment'};
        name='Add record?';
        numlines=1;
        defaultanswer={save_location,'None'};

        answer=inputdlg(prompt,name,numlines,defaultanswer);
        options.Resize='on';
        options.WindowStyle='normal';
        options.Interpreter='tex';

        if ~isempty(answer)
        tmp_struct=cell2struct(answer(2),'tmp',1); %save
        track.comment=(tmp_struct.tmp);
        tmp_struct=cell2struct(answer(1),'tmp',1); %save
        save_location=(tmp_struct.tmp)
            
        [xt yt]=find(figure_table==h);
        track.file_name=file_name(xt(end)).curr;
        add_record=[add_record ; track];
        eval(['save ' save_location ' add_record;']);
        end
    end

    
    return


%=========Open_song callback===========


function []=Open_song(hObject,eventdata,h,h_t,position_index)

global file_name figure_table fill_box save_location

    
    [songname, song_pathname, filterindex] = uigetfile( ...
       {'*.wav','WAV-files (*.wav)'; ...
        '*.*',  'All Files (*.*)'}, ...
        'Pick saved wav file', ...
        'MultiSelect', 'on');
    
    fill_box=0;
    if songname~=0
        if isempty(h)
            delete(h_t)
        end
    if isempty(h)
        file_name=[];
        file_name.curr=[song_pathname songname];
    else       
        xt=size(figure_table,1)+1;
        file_name(xt(end)).curr=[song_pathname songname];
    end
        [data,fs]=wavread([song_pathname songname]); 
        if nargin==4
            Call_spectrum_derivative(data,fs);
        else
            if~isempty(position_index)
                Call_spectrum_derivative(data,fs,position_index,[h,h_t]);
            else
                Call_spectrum_derivative(data,fs);
            end
        end
    end
    
save_location=[];
return

%=========Parameters callback===========


function []=Change_Parameters(hObject,eventdata,fs)
global param

load('Parameters');

prompt={'FFT data window','Advance FFT window by','spectrum range','Min frequency range of wiener entropy & amplitude','Max frequency range of wiener entropy & amplitude','Upper pitch bound','Lower pitch bound','Harmonic pitch must be lower than','Goodness of pitch is higher than','winer entropy is lower than','Pitch averaging by its goodness - 1 , else 0','Plot axes - off 0 or on 1?','Detection of silence low amplitude threshold?' };
name='Parameters';
numlines=1;
%defaultanswer={'9.27','1','11020','860','8600','11008','449','1800','100','-2','y','on','change'};
if param.pitch_averaging
    param.pitch_averaging='y';
else 
    param.pitch_averaging='n';
end

defaultanswer={num2str(param.window/44.1),num2str(param.winstep/44.1),num2str(param.spectrum_range*fs/param.window/2),num2str(param.min_freq_winer_ampl/param.pad*fs),num2str(param.max_freq_winer_ampl/param.pad*fs),num2str(param.up_pitch),num2str(param.low_pitch),num2str(param.pitch_HoP),num2str(param.gdn_HoP),num2str(param.up_wiener),num2str(param.pitch_averaging),num2str(param.initial_axes),num2str(param.low_amplitude_th),'Off','change'};

answer=inputdlg(prompt,name,numlines,defaultanswer);
options.Resize='on';
options.WindowStyle='normal';
options.Interpreter='tex';


if ~isempty(answer)
            
tmp_struct=cell2struct(answer(1),'tmp',1); %save
param.window=round(str2num(tmp_struct.tmp)*44.1);
tmp_struct=cell2struct(answer(2),'tmp',1); %save
param.winstep=round(str2num(tmp_struct.tmp)*44.1);
tmp_struct=cell2struct(answer(3),'tmp',1); %save
param.spectrum_range=round(str2num(tmp_struct.tmp)/fs*param.window*2);
tmp_struct=cell2struct(answer(4),'tmp',1); %save
param.min_freq_winer_ampl=round(str2num(tmp_struct.tmp)*param.pad/fs);
tmp_struct=cell2struct(answer(5),'tmp',1); %save
param.max_freq_winer_ampl=round(str2num(tmp_struct.tmp)*param.pad/fs);;
tmp_struct=cell2struct(answer(6),'tmp',1); %save
param.up_pitch=str2num(tmp_struct.tmp);    %?????????????????????????????????????????????????
tmp_struct=cell2struct(answer(7),'tmp',1); %save
param.low_pitch=str2num(tmp_struct.tmp);  %?????????????????????????????????????????????????
tmp_struct=cell2struct(answer(8),'tmp',1); %save
param.up_HoP=str2num(tmp_struct.tmp);
tmp_struct=cell2struct(answer(9),'tmp',1); %save
param.low_HoP=str2num(tmp_struct.tmp);
tmp_struct=cell2struct(answer(10),'tmp',1); %save
param.up_wiener=str2num(tmp_struct.tmp);
tmp_struct=cell2struct(answer(11),'tmp',1); %save
param.pitch_averaging=(tmp_struct.tmp(1));
if param.pitch_averaging=='y' | param.pitch_averaging=='Y'
    param.pitch_averaging=1;
else 
    param.pitch_averaging=0;
end
tmp_struct=cell2struct(answer(12),'tmp',1); %save
param.initial_axes=str2num(tmp_struct.tmp);
tmp_struct=cell2struct(answer(13),'tmp',1); %save
param.low_amplitude_th=str2num(tmp_struct.tmp);
do=(tmp_struct.tmp);


directoryname = uigetdir('c:\', 'Chose the directory where the parameters would be saved. It should be either in the matlab path or in the directory of the program');
 
while directoryname~=0
t=cd;
   if length(t)==length(directoryname)
       if min(t==directoryname)
           Parameters(param,do,directoryname);
           
           return
       end
   end


s = matlabpath;
p = 1;

while true
   t = strtok(s(p:end), pathsep);
   if length(t)==length(directoryname)
       if min(t==directoryname)
           Parameters(param,do,directoryname);
           
           return
       end
   end

   disp(sprintf('%s', t))
   p = p + length(t) + 1;
   

   if isempty(strfind(s(p:end),';')) break, end;
end



directoryname = uigetdir('c:\', 'Chose the directory where the parameters would be saved. It should be either in the matlab path or in the directory of the program');
end 

end
return


%=========Pan_x callback===========


function []=Pan_x(hObject,eventdata,h,h_t,bgcolor,m_spec_deriv,sound_bird, x_length , y_length , fs , m_AM, m_FM ,m_Entropy , m_Amplitude ,m_Freq,  m_Pitch_gdns , m_Pitch )

persistent pan_state

figure(h)
%cdata_index=round(adj_position([45,15]));
%cdata_color=ones(cdata_index(2),cdata_index(1))*0.8;
%cdata_color(:,:,1)=0.7;
%cdata_color(:,:,2)=0.8;
%cdata_color(:,:,3)=0.4;

ScreenSize=get(0,'ScreenSize');

if isempty(pan_state)
    pan_state=1;
    pan_sd(h,h_t,'xon')
    figure(h_t)
    hand=imread('hand_g.jpg');
    hand = imresiz(hand,[size(hand,1)*ScreenSize(4)/768 size(hand,2)*ScreenSize(3)/1024]);

    uicontrol('Style','pushbutton',...
        'Position',adj_position([ 3 163 15 20]),...   %18 140 45 15]),...
        'Interruptible','off',...
        'CData',hand,...
        'TooltipString','Pan',...
        'BusyAction','cancel',...
        'BackgroundColor',bgcolor,...
        'Callback',{@Pan_x,h,h_t,bgcolor,m_spec_deriv,sound_bird, x_length , y_length , fs , m_AM, m_FM ,m_Entropy , m_Amplitude ,m_Freq,   m_Pitch_gdns , m_Pitch });

%    set(h,'WindowButtonDownFcn',{@pan_fast_pan,h,h_t,size(m_spec_deriv,2)})
    
elseif pan_state
    pan_state=0;
    pan_sd(h,h_t,'off')
    figure(h_t)
    hand=imread('hand.jpg');
    hand = imresiz(hand,[size(hand,1)*ScreenSize(4)/768 size(hand,2)*ScreenSize(3)/1024]);

    uicontrol('Style','pushbutton',...
        'Position',adj_position([ 3 163 15 20]),...   %18 140 45 15]),...
        'Interruptible','off',...
        'CData',hand,...
        'TooltipString','Pan',...
        'BusyAction','cancel',...
        'BackgroundColor',bgcolor,...
        'Callback',{@Pan_x,h,h_t,bgcolor,m_spec_deriv,sound_bird, x_length , y_length , fs , m_AM, m_FM ,m_Entropy , m_Amplitude ,m_Freq,   m_Pitch_gdns , m_Pitch });
 %   set(h,    'WindowButtonDownFcn',{@get_pnts,m_spec_deriv,sound_bird, x_length , y_length ,h , h_t , fs , m_AM, m_FM ,m_Entropy , m_Amplitude ,m_Freq,   m_Pitch_gdns , m_Pitch   });
    

else
    hand=imread('hand_g.jpg');
    hand = imresiz(hand,[size(hand,1)*ScreenSize(4)/768 size(hand,2)*ScreenSize(3)/1024]);

    pan_state=1;
    pan_sd(h,h_t,'xon')
    figure(h_t)
    uicontrol('Style','pushbutton',...
     'Position',adj_position([ 3 163 15 20]),...   %18 140 45 15]),...
        'Interruptible','off',...
        'BusyAction','cancel',...
        'TooltipString','Pan',...
        'CData',hand,...
        'BackgroundColor',bgcolor,...
        'Callback',{@Pan_x,h,h_t,bgcolor,m_spec_deriv,sound_bird, x_length , y_length , fs , m_AM, m_FM ,m_Entropy , m_Amplitude ,m_Freq,   m_Pitch_gdns , m_Pitch });
  %  set(h,'WindowButtonDownFcn',{@pan_fast_pan,h,h_t,size(m_spec_deriv,2)})

end


return

%=========plot_function callback===========


function []=plot_feature(hObject,eventdata,x_length,y_length,plot_function,features_idn,index,h,h_t,cdata_color)

global bgcolor
load('Parameters');
ScreenSize=get(0,'ScreenSize');
 
          figure(h);
            h_line = findobj(gca,'Type','line');
            if ~isempty(h_line)
                
            for k=1:length(h_line)
                cdata_line=get(h_line(k),'Color');
                if min(cdata_line==features_idn(index).color)
                    delete(h_line(k));
                              figure(h_t);

                          uicontrol(gcf,'Style','pushbutton',...
                            'String', features_idn(index).str,...
                            'Interruptible','off',...
                            'CData',cdata_color,...
                            'BusyAction','cancel',...
                            'Callback',{@plot_feature,x_length,y_length,plot_function,features_idn,index,h,h_t,cdata_color },...
                            'Position',adj_position([15 features_idn(index).position*768/ScreenSize(4) 50 18]),...
                            'BackgroundColor',bgcolor);
                        return
                end
            end
            end
            



                          figure(h);
                          curr_axis=axis;
               %           m_const=(x_length)/0.732/(ScreenSize(3));
                          m_const=1;
                            modified_x_index=(1:length(plot_function))'*m_const;
%                          take_off=find(diff(round(modified_x_index))==0);
%                          modified_x_index(take_off)=[];
%                          plot_function(take_off)=[];
                          hold on
                          
                          if features_idn(index).position==20/768*ScreenSize(4); %M-Freq
                            load('Parameters');
                            tmp_plot_function=plot_function/1000;
                              plot((modified_x_index), tmp_plot_function  , 'Color',features_idn(index).color)
                          elseif features_idn(index).position==120/768*ScreenSize(4);
                            load('Parameters');
                            tmp_plot_function=plot_function/1000;
                              plot((modified_x_index), tmp_plot_function  ,'.', 'Color',features_idn(index).color)
                          
                          elseif features_idn(index).position==103/768*ScreenSize(4); %fm
                            load('Parameters');
                              tmp_plot_function=plot_function/90*y_length*param.fs/param.pad/1000;
                              plot((modified_x_index), tmp_plot_function  ,'.', 'Color',features_idn(index).color)
                          elseif features_idn(index).position==83/768*ScreenSize(4); %entropy
                            load('Parameters');
                            tmp_plot_function=plot_function+10;
                              tmp_plot_function=tmp_plot_function/10*y_length*param.fs/param.pad/1000;
                              plot((modified_x_index), tmp_plot_function  ,'.', 'Color',features_idn(index).color)
                          elseif features_idn(index).position==40/768*ScreenSize(4); %amplitude
                              
                             tmp_plot_function=plot_function-max(prctile_b(plot_function,5),0);
                              tmp_plot_function=tmp_plot_function/max(tmp_plot_function)*y_length*param.fs/param.pad/1000;
                              plot((modified_x_index), tmp_plot_function  , 'Color',features_idn(index).color)
                           
                          else %normelize
                              tmp_plot_function=plot_function-min(plot_function);
                              tmp_plot_function=tmp_plot_function/max(tmp_plot_function)*y_length*param.fs/param.pad/1000;
                              plot((modified_x_index), tmp_plot_function  , 'Color',features_idn(index).color)
                          end

                          figure(h_t);
                          tmp_cdata_color=cdata_color;
                          tmp_cdata_color(:,:,1)=features_idn(index).color(1);
                          tmp_cdata_color(:,:,2)=features_idn(index).color(2);
                          tmp_cdata_color(:,:,3)=features_idn(index).color(3);
                          
                          

                                uicontrol(gcf,'Style','pushbutton',...
                                'String', features_idn(index).str,...
                                'Interruptible','off',...
                                'BusyAction','cancel',...
                                'Callback',{@plot_feature,x_length,y_length,plot_function,features_idn,index,h,h_t,cdata_color },...
                                'Position',adj_position([15 features_idn(index).position*768/ScreenSize(4) 50 18]),...
                                'CData',tmp_cdata_color,...
                                'BackgroundColor',features_idn(index).color);

       
            
return

%=========Attach callback===========


function []=Attach(hObject,eventdata,h,h_t)

position_index=getappdata(h_t,'position_index');
bgcolor=[0.8 0.8 1];
position_h=[
          0 0.65 0.74 0.3 ; 
          0 0.4 0.74 0.3 ;
          0 0.1 0.74 0.3 ;];
position_h_t=[
          0.745 0.65 0.255 0.3 ; 
          0.745 0.4 0.255 0.3 ;
          0.745 0.1 0.255 0.3 ;];
ScreenSize=get(0,'ScreenSize');
ScreenSize=[ScreenSize(3:4) ScreenSize(3:4)];
position_h=position_h.*(ones(3,1)*ScreenSize);
position_h_t=position_h_t.*(ones(3,1)*ScreenSize);
%position_h_t(:,3)=0.255*1024;
position_h_t(:,1)=ScreenSize(3)-position_h_t(1,3);
%position_h(:,3)=position_h_t(:,1)-0.005*1024;
position_h(:,3)=position_h_t(:,1)-0.005*ScreenSize(3);

    set(h,'Units','pixel',...
    'position', position_h(position_index- floor((position_index)/3)*3+1,:));

    set(h_t,'Units','pixel',...
    'position', position_h_t(position_index- floor((position_index)/3)*3+1,:));
    figure(h_t)
    figure(h)
return


%=========plot axis callback===========


function []=plot_axis(hObject,eventdata,h,h_t,m_spec_deriv)

global cutoff   axis_x1 axis_x2 axis_y1 axis_y2 bgcolor
figure(h)
 v=getappdata(h, 'axes_index');
load('Parameters');
            h_line = findobj(gca,'Type','line');
            if ~isempty(h_line)
                
            for k=1:length(h_line)
                cdata_line=get(h_line(k),'Color');
            end
            end

onoff=get(v,'Visible');
if onoff(2)=='f'
    
        figure(h);
        old_axis=axis;

        h_line = findobj(gca,'Type','line');
       if ~isempty(h_line) 
        for j=1:length(h_line)
            save_line(j)=get(h_line(j));

        end
       end
       
        h_obj = findobj(h);

        delete(h_obj(2:end))
        
        v=axes('position', [0.03 0.11 0.97 0.89]);
        setappdata(h, 'axes_index',v);

% axes('position', [0 0 1 1])
        imagesc((1:size(m_spec_deriv,2)),param.fs*(1:size(m_spec_deriv,1))/1000/param.pad,m_spec_deriv);     
        %axis([axis_x1 axis_x2 param.fs*axis_y1/param.pad/1000 param.fs*axis_y2/param.pad/1000]); 
        axis(old_axis);
        colormap(bone); set(gca,'YDir','normal'); caxis([-cutoff cutoff])
        axis on
        ylabel('KHz','Units','normalized','Position',[0,0.05]);
        xlabel('ms','Units','normalized','Position',[0.05,0]);
       if ~isempty(h_line) 
        for j=1:length(h_line)
            hold on

            plot(save_line(j).XData,save_line(j).YData,'Color',save_line(j).Color,'Marker',save_line(j).Marker)
        end
       end
 
else
        figure(h);
        old_axis=axis;
        h_line = findobj(gca,'Type','line');
       if ~isempty(h_line) 
        for j=1:length(h_line)
            save_line(j)=get(h_line(j));

        end
       end
        h_obj = findobj(h);
        delete(h_obj(2:end))

        v=axes('position', [0 0 1 1]);
        setappdata(h, 'axes_index',v);
        imagesc((1:size(m_spec_deriv,2)),param.fs*(1:size(m_spec_deriv,1))/1000/param.pad,m_spec_deriv);     
        axis(old_axis);
        %axis([axis_x1 axis_x2 param.fs*axis_y1/param.pad/1000 param.fs*axis_y2/param.pad/1000]); axis off;
        axis off
        colormap(bone); set(gca,'YDir','normal'); caxis([-cutoff cutoff])
        
        if ~isempty(h_line) 
        for j=1:length(h_line)
            hold on
            plot(save_line(j).XData,save_line(j).YData,save_line(j).Marker,'Color',save_line(j).Color,save_line(j).Marker)
        end
        end


end
return




%=========play sound callback===========


function []=play_sound_bird(hObject,eventdata,x,sound_bird,fs,ration_sd,h,h_t)

           figure(h);
            play_all=1;
           h_line = findobj(gca,'Type','line');
            for k=1:length(h_line)
                xdata=get(h_line(k),'XData');
                if length(xdata)==2
                    play_all=0;
                end
            end

if play_all
    wavplay(sound_bird,fs);
else
    wavplay(sound_bird(floor(x(1)*ration_sd):ceil(x(2)*ration_sd)),fs)
end
return



%=========Rhythm callback===========


function []=rhythm_spectral_derivatives(hObject,eventdata,fs,  m_AM, m_FM ,m_Entropy , m_Amplitude , m_Freq , m_Pitch_gdns , m_Pitch )

   prompt={'Amplitude','AM','FM','Entropy','Mean Frequency','Pitch','Pitch goodness'};
   name='Rhythm on';
   numlines=1;
   defaultanswer={'y','y','y','y','y','y','y'};

   options.Resize='on';
   options.WindowStyle='normal';
   options.Interpreter='tex';
 
   answer=inputdlg(prompt,name,numlines,defaultanswer);
 
   tmp_struct=cell2struct(answer(1),'tmp',1); %save
   do_amplitude=(tmp_struct(1).tmp);
   
   tmp_struct=cell2struct(answer(2),'tmp',1); %save
   do_AM=(tmp_struct(1).tmp);

   tmp_struct=cell2struct(answer(3),'tmp',1); %save
   do_FM=(tmp_struct(1).tmp);
 
   tmp_struct=cell2struct(answer(4),'tmp',1); %save
   do_entropy=(tmp_struct(1).tmp);

   tmp_struct=cell2struct(answer(5),'tmp',1); %save
   do_freq=(tmp_struct(1).tmp);

   tmp_struct=cell2struct(answer(6),'tmp',1); %save
   do_pitch_good=(tmp_struct(1).tmp);

   tmp_struct=cell2struct(answer(7),'tmp',1); %save
   do_pitch=(tmp_struct(1).tmp);

   %web('sam.htm')
param.pad=8192*4;
fft_length=param.pad;
f = (0:(fft_length/2))/fft_length;
f_max=60;
index_f_max=find((f*1000)>f_max);
index_f_max=index_f_max(1);
f=f(1:index_f_max)*1000;%f times sampling frequancy
S=zeros(size(f));
index=0;
title_is=['Rhythm of: '];
if do_amplitude(1)=='y'
    data=m_Amplitude;
    [E,v]=dpss(length(data),1.5,1);
    J1=(fft(data.*(E(:,1)./max(E(:,1))),param.pad))/(param.pad);
    S_tmp=J1.*conj(J1); S=S_tmp(floor(end/2):end);
    index=index+1;
    title_is=[title_is 'Amplitude,'];
end
if do_AM(1)=='y'
    data=m_AM;
    [E,v]=dpss(length(data),1.5,1);
    J1=(fft(data.*(E(:,1)./max(E(:,1))),param.pad))/(param.pad);
    S_tmp=J1.*conj(J1); S=S_tmp(floor(end/2):end);
    index=index+1;
    title_is=[title_is 'AM,'];
end
if do_FM(1)=='y'
    data=m_FM;
    [E,v]=dpss(length(data),1.5,1);
    J1=(fft(data.*(E(:,1)./max(E(:,1))),param.pad))/(param.pad);
    S_tmp=J1.*conj(J1); S=S_tmp(floor(end/2):end);
    index=index+1;
    title_is=[title_is 'FM,'];
end
if do_entropy(1)=='y'
    data=m_Entropy;
    [E,v]=dpss(length(data),1.5,1);
    J1=(fft(data.*(E(:,1)./max(E(:,1))),param.pad))/(param.pad);
    S_tmp=J1.*conj(J1); S=S_tmp(floor(end/2):end);
    index=index+1;
    title_is=[title_is 'Entropy,'];
end
if do_freq(1)=='y'
    data=m_Freq;
    [E,v]=dpss(length(data),1.5,1);
    J1=(fft(data.*(E(:,1)./max(E(:,1))),param.pad))/(param.pad);
    S_tmp=J1.*conj(J1); S=S_tmp(floor(end/2):end);
    index=index+1;
    title_is=[title_is 'Mean Frequency,'];
end
if do_pitch_good(1)=='y'
    data=m_Pitch_gdns;
    [E,v]=dpss(length(data),1.5,1);
    J1=(fft(data.*(E(:,1)./max(E(:,1))),param.pad))/(param.pad);
    S_tmp=J1.*conj(J1); S=S_tmp(floor(end/2):end);
    index=index+1;
    title_is=[title_is 'Pitch goodness,'];
end
if do_pitch(1)=='y'
    data=m_Pitch;
    [E,v]=dpss(length(data),1.5,1);
    J1=(fft(data.*(E(:,1)./max(E(:,1))),param.pad))/(param.pad);
    S_tmp=J1.*conj(J1); S=S_tmp(floor(end/2):end);
    index=index+1;
    title_is=[title_is 'Pitch,'];
end
if index~=0
    S=S/index;
    figure
    plot(f,S(1:length(f)))
    xlabel('Hz');
    title(title_is(1:end-1))

else
    errordlg('No feature was chosen')
end

return

%=========asosicate paning callback===========


function []=pan_fast_pan(hObject,eventdata,h,h_t,data_length)

%ScreenSize=get(0,'ScreenSize');
%curr_position=get(gcf,'Position');

%A=selectmoveresize
curr_point=(get(h_t,'XLim'))
curr_point=mean(curr_point); %get(A.Handles,'XLim'));


figure(h_t)

h_line = findobj(gca,'Type','line');
delete(h_line);
h_l=plot([curr_point(1) curr_point(1)],[0 500],'r');


return

