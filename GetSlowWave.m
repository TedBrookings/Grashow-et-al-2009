function [Freq, varargout] = GetSlowWave(DeltaT, V0, V1, PlotSubject)
%[Freq, NumSigma, SlowC] = GetSlowWave(DeltaT, V0, V1, PlotSubject)
% calculates the slow wave frequency and sigma
%   INPUT PARAMETERS:
%    -DeltaT is time step (s)
%    -V0 is first voltage trace (arbitrary units)
%    OPTIONAL:
%     -V1 is second voltage trace (arbitrary units)
%     -PlotSubject  set to 1/true to plot.  Defaults to false
%   OUTPUT PARAMETERS:
%     -Freq is dominant slow wave frequency in Hz
%      OPTIONAL:
%      -NumSigma is z-score of slow wave frequency correlation,
%       compared to average
%      -SlowC is autocorrelation at slow wave period

%May want to build in ability to use pmtm or coh_mt.m

CCutoff = 0.05;
%First calculate the autocorrelations:
NumAutoCorr = round(length(V0) / 3);  %need to see at least 3 waves!
C = xcorr(zscore(V0), NumAutoCorr, 'unbiased');
C = C((NumAutoCorr+1):end);
f = (1.0 / DeltaT) ./ (1:length(C));

if(nargin >= 3)
  if(length(V1) == length(V0))
    C1 = xcorr(zscore(V1), NumAutoCorr, 'unbiased');
    C1 = C1((NumAutoCorr+1):end);
    C = 0.5 * (C + C1);
    if(nargin == 3)
      PlotSubject = false;
    end
  elseif(nargin == 3)
    PlotSubject = V1;
  else
    error('Invalid value for V1')
  end
else
  PlotSubject = false;
end

%Find the slow wave peak.
[MaxInd, MaxC, IndStart, IndStop] = FindPeak(C);
Freq = f(MaxInd);
if(length(Freq) == 0)
  Freq = 0;
  MaxC = 0;
  IndStart = 1;
  IndStop = 1;
elseif(MaxC < CCutoff)
  TempStart = IndStart;
  IndStart = find(C < CCutoff, 1);
  IndStart = find(C((IndStart + 1):end) >= CCutoff, 1) + IndStart;
  if(length(IndStart) == 0)
    IndStart = TempStart;
  else
    IndStop = find(C((IndStart+1):end) < CCutoff, 1) + IndStart;
    if(length(IndStop) == 0)
      IndStop = length(C);
    end
    [MaxC, MaxInd] = max(C(IndStart:IndStop));
    MaxInd = MaxInd + IndStart - 1;
    Freq = f(MaxInd);
  end
end
if(DoPlot(PlotSubject))
  if(ischar(PlotSubject) & length(PlotSubject) > 0)
    TitleStr = ['Slow-wave Correlation for ', PlotSubject];
  else
    TitleStr = 'Slow-wave Correlation';
  end
  h = NamedFigure(TitleStr);
  set(h, 'WindowStyle', 'docked');
  hold off
  plot(f, C, 'b.')
  if(Freq > 0)
    hold on
    CRange = [min(C), MaxC];
    plot([f(IndStart), f(IndStart)], CRange, 'g-')
    plot([f(IndStop), f(IndStop)], CRange, 'g-')
    plot(Freq, MaxC, 'ro', 'MarkerFaceColor', 'r')
    hold off
    fStop = max([2 * Freq, f(IndStart) * 1.1]);
    xlim([0, fStop])
  else
    xlim([0, 5])
  end
  xlabel('Frequency (Hz)', 'FontSize', 18);
  ylabel('AutoCorrelation', 'FontSize', 18);
  title(RealUnderscores(TitleStr), 'FontSize', 18);
end

if(MaxC < CCutoff)
  Freq = 0;
end

%Calculate NumSigma if needed
if(nargout > 1)
  if(MaxC <= CCutoff)
    NumSigma = 0;
  else
    NumSigma = abs(MaxC - mean(C)) / std(C);
  end
  if(nargout == 2)
    varargout = {NumSigma};
  else
    varargout = {NumSigma, MaxC};
  end
else
  varargout = {};
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotVar = DoPlot(PlotSubject)
if(ischar(PlotSubject))
  PlotVar = true;
else
  PlotVar = PlotSubject;
end
return
