function C = corr_rate_smooth(spktrain, trange, twidth, tres, varargin)
% CORR_RATE_SMOOTH
% 
% Usage:  C = corr_rate_smooth(sptrain, trange, twidth, tres, varargin)
% 
% sptrain   a cell array of spike trains. Each spike train is a an
%           array of spike times in msec. There should be at least
%           two spike trains. If there are N trains, all possible
%           pairwise distances are returned in an NxN upper
%           triangular matrix, D.
%
% trange    the time range in msec to consider in the spike trains.  
%
% twidth    The duration corresponding to 1 std dev. of the
%           Gaussian filter used for smoothing the spike trains.
%
% tres      the time resolution in msec to keep for doing the 
%           convolution.
%

if nargin > 4
  do_shuffle = varargin{1};
  % Reset random number generator state so we get the same results
  % each time.
  rand('state',137)

else
  do_shuffle = 0;
end

n_trains = length(spktrain);
C = zeros(n_trains,n_trains);

trial_dur = trange(2) - trange(1);
n_points = 1 + round(trial_dur/tres);

if do_shuffle
  % Shift (circularly) each trial a random amount to generate a shuffled estimate.

  shift_min = fix(100.0/tres); % 100ms min shift.
  shift_max = fix(min(500.0,trial_dur-100.0)/tres);
  if shift_max <= shift_min | trial_dur <= 100.0
    shift_min = fix(trial_dur/(4.0*tres));
    shift_max = fix(trial_dur/(2.0*tres));
  end
  dshift = shift_max - shift_min;
  for i=1:n_trains
    shift_rand(i) = fix(shift_min + rand(1,1)*dshift);
    posneg = rand(1,1);
    if posneg < 0.5
      posneg = -1.0;
    else
      posneg = 1.0;
    end
    shift_rand(i) = posneg*shift_rand(i);
  end
end

% Number of multiples of 1 std dev. at which to truncate Gaussian window.
trunc_gauss = 4;
width = twidth/tres;
n_filt_pts = 1 + round(2*trunc_gauss*width);
%exp_filter = exp(-(0:n_filt_pts-1)*r);
gauss_filter = gausswin(n_filt_pts, 1/width);

figure;
hold on;

for i=1:n_trains
  sptrain = (spktrain{i});
  sptimes1_idx = find(sptrain <= trange(2) & sptrain >= trange(1));
  n_spikes(i) = length(sptimes1_idx);
  sptimes1 = sptrain(sptimes1_idx);
  spseries1_idx = 1 + round((sptimes1 - trange(1))/tres);
  spseries1 = zeros(1,n_points);
  spseries1(spseries1_idx) = 1;    
  spseries1 = conv(gauss_filter,spseries1);
  plot(spseries1);
  if do_shuffle
    spseries1 = circshift(spseries1',shift_rand(i))';
  end
  sprate1_mean = mean(spseries1);

  for j=i+1:n_trains
    sptrain = (spktrain{j});
    % Make a binary time series of resolution tres out of discrete spike times.
    sptimes2_idx = find(sptrain <= trange(2) & sptrain >= trange(1));
    sptimes2 = sptrain(sptimes2_idx);
    spseries2_idx = 1 + round((sptimes2 - trange(1))/tres);
    spseries2 = zeros(1,n_points);
    spseries2(spseries2_idx) = 1;    
    spseries2 = conv(gauss_filter,spseries2);
    if do_shuffle
      spseries2 = circshift(spseries2',shift_rand(j))';
    end
    sprate2_mean = mean(spseries2);
    
    denom = sqrt(sum((spseries1 - sprate1_mean).^2)*sum((spseries2 - sprate2_mean).^2));
    if denom ~= 0
      C(i,j) = sum((spseries1 - sprate1_mean).*(spseries2 - sprate2_mean))/denom;
    else
      C(i,j) = 0.0;
    end
    
  end
end


