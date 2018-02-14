function Ramp = GetRampData(varargin)
% Ramp = GetRampData(FileName, CellName, PlotVar)
%        ---OR---
% Ramp = GetRampData(t, v, I, PlotVar)
% if no parameters are passed, a structure filled with NaNs is returned.
%  INPUT PARAMETERS:
%   -t is time in ms
%   -v is voltage in mV
%   -I is current in nA
%   -FileName is the name of an .abf file with IV data
%   -CellName is the name of the requested cell within the
%      file specified by FileName (e.g. "top" or "bot")
%    OPTIONAL:
%     -PlotVar set to 1/true to plot the result
%     (defaults to 0/false)
%  OUTPUT PARAMETERS:
%   -Ramp: structure with information related to current ramp data
%     -Ramp.Analyze:  structure produced by AnalyzeWaveform3.m
%     -Ramp.SpikeHeight:  voltage of spike peak, minus voltage at
%          point of maximum curvature before spike
%     -Ramp.SpikeWidth:  time of spike peak, minus time at point of
%          maximum curvature before spike
%     -Ramp.TenMs: structure with I/V data for point 10 ms before spike peak
%       -Ramp.TenMs.V:  voltage in mV
%       -Ramp.TenMs.I:  current in nA
%     There are several other structures, with hopefully
%      self-descriptive names.  They all contain, a .t, .V and .I
%      field, corresponding to the time that they are sampled, the
%      voltage, and the current.  They may also have a fourth field
%      (e.g. K for maximum curvature).  They are:
%        -Ramp.MaxV, Ramp.MaxDeriv, Ramp.MinDeriv,
%         Ramp.PreMaxCurve, Ramp.PostMaxCurve

Ramp.MaxV.V = NaN;
Ramp.MaxV.t = NaN;
Ramp.MaxV.I = NaN;
Ramp.MaxDeriv.V = NaN;
Ramp.MaxDeriv.DV = NaN;
Ramp.MaxDeriv.t = NaN;
Ramp.MinDeriv.V = NaN;
Ramp.MinDeriv.I = NaN;
Ramp.MinDeriv.DV = NaN;
Ramp.MinDeriv.t = NaN;
Ramp.MinDeriv.I = NaN;
Ramp.PreMaxCurve.t = NaN;
Ramp.PreMaxCurve.V = NaN;
Ramp.PreMaxCurve.K = NaN;
Ramp.PreMaxCurve.I = NaN;
Ramp.PostMaxCurve.t = NaN;
Ramp.PostMaxCurve.V = NaN;
Ramp.PostMaxCurve.K = NaN;
Ramp.PostMaxCurve.I = NaN;

if(nargin < 1)
  Okay = false;
  BadNumArg = false;
elseif(ischar(varargin{1}))   %specified FileName
  if(nargin < 2 | nargin > 3)
    BadNumArg = true;
  else
    BadNumArg = false;
    FileName = varargin{1};
    CellName = varargin{2};
    if(nargin == 3)
      PlotVar = varargin{3};
    else
      PlotVar = false;
    end
    [t, v, I, Okay] = LoadCleanIntrinsicData(FileName, CellName);
  end
else
  if(nargin < 3 | nargin > 4)
    BadNumArg = true;
  else
    BadNumArg = false;
    Okay = true;
    t = varargin{1};
    v = varargin{2};
    I = varargin{3};
    if(nargin == 4)
      PlotVar = varargin{4};
    else
      PlotVar = false;
    end
  end
end
if(BadNumArg)
  error('Incorrect number of input arguments.  Run "help GetRampData"')
elseif(~Okay)
  Ramp.SpikeHeight = NaN;
  Ramp.SpikeWidth = NaN;
  Ramp.TenMs.I = NaN;
  Ramp.TenMs.V = NaN;
  Ramp.Analyze = NaN;
  return
end

NoShape = false;
FirstOnly = true;

Analyze = AnalyzeWaveform3(t, v, PlotVar, NoShape, FirstOnly);
FNames = fieldnames(Ramp);
for n = 1:length(FNames)
  FName = FNames{n};
  FirstSpikeT = Analyze.Spike.(FName).t;
  if(length(FirstSpikeT) ~= 1)
    AnalyzeWaveform3(t, v, true, NoShape, FirstOnly);
    FirstSpikeT
    error('Incorrect number of spikes found')
  end
  [MinDiff, n] = min(abs(t - FirstSpikeT));
  Ramp.(FName).I = I(n);
  
  SubFNames = fieldnames(Ramp.(FName));
  for m = 1:length(SubFNames)
    SubFName = SubFNames{m};
    if(strcmp(SubFName, 'I'))
      continue
    elseif(strcmp(SubFName, 't'))
        Ramp.(FName).(SubFName) = Analyze.Spike.(FName).(SubFName);
    else
        Ramp.(FName).(SubFName) = Analyze.Spike.(FName).(SubFName).List;
    end
  end
end
Ramp.SpikeHeight = Ramp.MaxV.V - Ramp.PreMaxCurve.V;
Ramp.SpikeWidth = Analyze.Spike.PostMinV.t - Ramp.PreMaxCurve.t;

tTenMs = Ramp.MaxV.t - 10;
[MinVal, Ind] = min(t - tTenMs);
Ramp.TenMs.I = I(Ind);
Ramp.TenMs.V = v(Ind);

Ramp.Analyze = Analyze;
return