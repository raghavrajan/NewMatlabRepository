global SS_MAX_ELECTRODE_VOLTAGE;
global SS_AD_MAX_VALUE;	 
global SS_TICKS_PER_VOLT;

SS_MAX_ELECTRODE_VOLTAGE = 0.01; % in volts
SS_AD_MAX_VALUE = 32767;	 % for 16 bit samples
SS_TICKS_PER_VOLT = SS_AD_MAX_VALUE/SS_MAX_ELECTRODE_VOLTAGE;

% This converts volts to 16 bit AD sample values.
%multiplier = ticks_per_volt/gain_external; 
