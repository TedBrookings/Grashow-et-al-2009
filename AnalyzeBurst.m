function Burst = AnalyzeBurst(Spike, SlowWave, t)
%   -Burst:  structure with burst information
%      -Burst.Freq is mean burst frequency (Hz)
%      -Burst.SpikeFreq is mean within-burst spike frequency (Hz)
%      -Burst.DutyCycle is the average burst duration/period
%      -Burst.Times is a plain list of burst times (ms)
%      -Burst.Durations is a list of burst durations (ms)
%      -Burst.SpikesPerBurst is a list of spikes per burst
%      -Burst.SpikeFrequencies is a list of spike frequencies (Hz)
%      -Burst.InterBurstIntervals is a list of intervals between bursts (ms)
%      -Burst.StartInds is a list of indices where bursts begin
%      -Burst.StopInds is a list of indices where bursts end

if(Spike.Freq <= 0 | (SlowWave.NumOnlySpike >= SlowWave.NumBurst))
  %This is not a bursting waveform
  Burst.Times = [];
  Burst.Durations = [];
  Burst.SpikesPerBurst = [];
  Burst.SpikeFrequencies = [];
  Burst.InterBurstIntervals = [];
  Burst.StartInds = [];
  Burst.StopInds = [];
  Burst.Freq = 0;
  Burst.SpikeFreq = 0;
  Burst.DutyCycle = 0;
  return
end

BurstInds = SlowWave.MinInds;
NumBurst = length(BurstInds) - 1;
CosPhase = cos(SlowWave.Phases);
t_Off = -Inf;

Burst.Times = zeros(NumBurst, 1);
Burst.Durations = zeros(NumBurst, 1);
Burst.SpikesPerBurst = zeros(NumBurst, 1);
Burst.SpikeFrequencies = zeros(NumBurst, 1);
Burst.InterBurstIntervals = zeros(NumBurst - 1, 1);
Burst.StartInds = zeros(NumBurst, 1);
Burst.StopInds = zeros(NumBurst, 1);

for BurstNum = 1:NumBurst
  Last_Off = t_Off;
  n1 = BurstInds(BurstNum);
  n2 = BurstInds(BurstNum + 1);
  n_On = find(CosPhase(n1:n2) < 0, 1) + n1 - 1;
  n_Off = find(CosPhase(n_On:n2) < 0, 1, 'last') + n_On - 1;

  SpikeInds = find(Spike.n1List > n1 & Spike.n1List <= n2);
  Spike_n1s = Spike.n1List(SpikeInds);
  Spike_n2s = Spike.n2List(SpikeInds);

  if(length(SpikeInds) > 1 && ...
     (Spike_n1s(1) < n_On || Spike_n2s(end) > n_Off))
    %either expand burst to include spikes, or
    %  remove non-burst spikes from list
    InBurst = find(Spike_n2s >= n_On & Spike_n1s <= n_Off);
    if(length(InBurst) > 1)
      if(Spike_n1s(InBurst(1)) < n_On)
	n_On = Spike_n1s(InBurst(1));
      end
      if(Spike_n2s(InBurst(end)) > n_Off)
	n_Off = Spike_n2s(InBurst(end));
      end
      MaxDiff = max(diff(Spike_n1s(InBurst)));
      while(InBurst(1) > 1)
	NewDiff = Spike_n1s(InBurst(1)) - Spike_n2s(InBurst(1) - 1);
	if(NewDiff < MaxDiff)
	  InBurst = [InBurst(1)-1; InBurst(:)];
	  n_On = Spike_n1s(InBurst(1));
	  NewDiff = Spike_n1s(InBurst(2)) - n_On;
	  if(NewDiff > MaxDiff)
	    MaxDiff = NewDiff;
	  end
	else
	  break
	end
      end
      NumOrig = length(SpikeInds);
      while(InBurst(end) < NumOrig)
	NewDiff = Spike_n1s(InBurst(end)+1) - Spike_n2s(InBurst(end));
	if(NewDiff < MaxDiff)
	  InBurst = [InBurst(:); InBurst(end)+1];
	  n_Off = Spike_n2s(InBurst(end));
	  NewDiff = Spike_n1s(InBurst(end)) - Spike_n1s(InBurst(end-1));
	  if(NewDiff > MaxDiff)
	    MaxDiff = NewDiff;
	  end
	else
	  break
	end
      end      
    end
    SpikeInds = SpikeInds(InBurst);
  end

  NumSpikes = length(SpikeInds);
  Burst.StartInds(BurstNum) = n_On;
  Burst.StopInds(BurstNum) = n_Off;
  t_On = t(n_On);
  t_Off = t(n_Off);
  Burst.Times(BurstNum) = t_On;
  Burst.Durations(BurstNum) = t_Off - t_On;

  Burst.SpikesPerBurst(BurstNum) = NumSpikes;
  if(NumSpikes > 2)
    Interval = .001 * (Spike.Times(SpikeInds(end)) ...
		       - Spike.Times(SpikeInds(1)));  %sec
    Burst.SpikeFrequencies(BurstNum) = (NumSpikes-1) / Interval;
  else
    Burst.SpikeFrequencies(BurstNum) = NaN;
  end
  if(BurstNum > 1)
    Burst.InterBurstIntervals(BurstNum) = Last_Off - t_On;
  end
end

if(NumBurst >= 3)
  Periods = .001 * (Burst.Times(2:end) - Burst.Times(1:(end-1)));
  Burst.Freq = mean(1.0 ./ Periods);
else
  Burst.Freq = 0;
end
Ind = find(isfinite(Burst.SpikeFrequencies));
Burst.SpikeFreq = mean(Burst.SpikeFrequencies(Ind));
Burst.DutyCycle = sum(Burst.Durations) / (t(BurstInds(end)) - t(BurstInds(1)));
return