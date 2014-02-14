function [add_record]=taper_spectrum_derivative(sound,position_index,cutoff, full_view , fs)

%DERIV   [S, S_f, S_t]=DERIV(TS,NW,K,PAD,WINDOW,WINSTEP)
%	 uses the derivative estimates to calculate the spectrum's frequency
%	 and time derivatives
% 
%        S: estimated spectrum; S_f: estimated frequency derivative; 
%        S_t: estimated time derivative
%        NW: time bandwidth parameter (e.g. 3)
%        K : number of tapers kept, approx. 2*NW-1
%        pad: length to which data will be padded (preferably power of 2
%        window: time window size
%        winstep: distance between centers of adjascent time windows
%        position_index: position on the screen gets values 1-3
%        cutoff - trunct image
%        full_view 1 for no axis 2 for smaller image with axis
%        
%        Writen by Sigal Saar August 08 2005

global cutoff quit_ntn cutoff_value


if nargin<1
    error('No sound file');
end

if nargin<2
    position_index=3;
end
if nargin<3
    cutoff=50*758400000;
end
if nargin<4
    full_view=2;
end

if nargin<5
    fs=44100;
end
    pad=1024;
    window=409;
    winstep=44;
cutoff_value=0.5;
[m_spec_deriv , m_AM, m_FM ,m_Entropy , m_amplitude , m_Freq, m_PitchGoodness , m_Pitch]=deriv(sound,pad,window,winstep,fs);

trunk(m_spec_deriv,  position_index , full_view ,[]  , m_AM, m_FM ,m_Entropy , m_amplitude ,m_Freq, m_PitchGoodness , m_Pitch , fs);

    