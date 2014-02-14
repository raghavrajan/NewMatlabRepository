function C = corr_isi_rate(spktrain, trange, tres, varargin)
% CORR_ISI_RATE
% 
% Usage:  C = corr_isi_rate(sptrain, trange, tres)
% 
% sptrain   a cell array of spike trains. Each spike train is a an
%           array of spike times in msec. There should be at least
%           two spike trains. If there are N trains, all possible
%           pairwise distances are returned in an NxN upper
%           triangular matrix, D.
%
% trange    the time range in msec to consider in the spike trains.  
%
% tres      the time resolution in msec to keep for doing the 
%           convolution.
%

if nargin > 3
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

figure;
hold on;
for i=1:n_trains
  sptrain = (spktrain{i});
  sptimes1_idx = find(sptrain <= trange(2) & sptrain >= trange(1));
  sptimes1 = sptrain(sptimes1_idx);
  nspikes1 = length(sptimes1);
  sprate1 = 1.0./diff(sptimes1');
  spseries1_idx = 1 + round((sptimes1 - trange(1))/tres);
  spseries1 = zeros(1,n_points);
  for k=1:nspikes1-1
    spseries1(1,spseries1_idx(k):spseries1_idx(k+1)) = sprate1(k); % rate in kHz
  end
  if do_shuffle
    spseries1 = circshift(spseries1',shift_rand(i))';
  end
  sprate1_mean = mean(spseries1);
    
  if (mod(i,2) == 0)
      plot(spseries1,'b'); 
  else
      plot(spseries1,'r');
  end
  hold on;
%   title(['Instantaneous spike rate, Trial ', num2str(i)])
%   pause
  
  for j=(i+1):n_trains
    sptrain = (spktrain{j});
    % Make a binary time series of resolution tres out of discrete spike times.
    sptimes2_idx = find(sptrain <= trange(2) & sptrain >= trange(1));
    sptimes2 = sptrain(sptimes2_idx);
    nspikes2 = length(sptimes2);
    sprate2 = 1.0./diff(sptimes2');
    spseries2_idx = 1 + round((sptimes2 - trange(1))/tres);
    spseries2 = zeros(1,n_points);
    for l=1:nspikes2-1
      spseries2(1,spseries2_idx(l):spseries2_idx(l+1)) = sprate2(l); % rate in kHz
    end
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


