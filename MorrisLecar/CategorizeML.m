function [Cat, CatString] = CategorizeML(Analysis)
%[Cat, CatString] = CategorizeML(Analysis)
% categorizes results from AnalyzeML and puts into one of four
% mutually-exclusive categories:
%  0 - Silent
%  1 - Spiking
%  2 - Slow wave
%  3 - Slow-wave with Spikes

SlowWaveCutoff = 2.0;  %Minimum amplitude to be "significant"


SpikeFreq = Analysis.CellReal.Spike.Freq;
BurstFreq = Analysis.CellReal.Burst.Freq;
HalfCenter = Analysis.CellReal.HalfCenter;

RealSlowWaveAmp = mean(Analysis.CellReal.SlowWave.Amplitudes);
if(isnan(RealSlowWaveAmp) | RealSlowWaveAmp < SlowWaveCutoff)
  RealSlowWaveAmp = 0;
end
ModelSlowWaveAmp = mean(Analysis.ModelSlow.Amplitudes);
if(isnan(ModelSlowWaveAmp) | ModelSlowWaveAmp < SlowWaveCutoff)
  ModelSlowWaveAmp = 0;
end

if(SpikeFreq <= 0)
  if(ModelSlowWaveAmp == 0)
    Cat = 0;
    CatString = 'Silent';
  else
    Cat = 2;
    CatString = 'Model slow-wave, GM inhibited';
  end
elseif(ModelSlowWaveAmp > 0 && HalfCenter)
  Cat = 3;
  CatString = 'Bursting';
else
  Cat = 1;
  CatString = 'GM spiking, Model inhibited';
end

return