function varargout = FindPeak(Y, StartInd, Direction, DerivLen)
if(nargin < 4 || length(DerivLen) == 0)
  DerivLen = 10;
if(nargin < 3 || length(Direction) == 0)
  Direction = +1;
end
if(nargin < 2 || length(StartInd) == 0)
  if(Direction < 0)
    StartInd = length(Y);
  else
    StartInd = 1;
  end
end

if(Direction < 0)
  Deriv = Y((DerivLen+1):StartInd) - Y(1:(StartInd-DerivLen));
  Trivial = 0.2 * median(abs(Deriv))
  Neg1 = find(Deriv < -Trivial, 1, 'last');
  Pos1 = find(Deriv(1:Neg1) > Trivial, 1, 'last');
  Neg2 = find(Deriv(1:Pos1) < -Trivial, 1, 'last');
  Pos2 = find(Deriv(1:Neg2) > Trivial, 1, 'last') - DerivLen;
  LookStart = Pos2;
  LookStop = Pos1;
  [PeakVal, PeakInd] = max(Y(LookStart:LookStop));
  PeakInd = PeakInd + LookStart - 1;
else
  Deriv = Y((DerivLen + StartInd):end) - Y(StartInd:(end-DerivLen));
  Trivial = 0.2 * median(abs(Deriv));
  Neg1 = find(Deriv < -Trivial, 1);
  Pos1 = find(Deriv(Neg1:end) > Trivial, 1) + Neg1 - 1;
  Neg2 = find(Deriv(Pos1:end) < -Trivial, 1) + Pos1 - 1;
  Pos2 = find(Deriv(Neg2:end) > Trivial, 1) + Neg2 - 1;
  LookStart = Pos1;
  LookStop = Pos2;
  [PeakVal, PeakInd] = max(Y(LookStart:LookStop));
  PeakInd = PeakInd + LookStart + StartInd - 2;
end
%{
figure
plot(Deriv)
hold on
plot(Neg1, Deriv(Neg1), 'go')
plot(Pos1, Deriv(Pos1), 'ro')
plot(Neg2, Deriv(Neg2), 'gx')
plot(Pos2, Deriv(Pos2), 'rx')
hold off
%}
varargout = {PeakInd, PeakVal, LookStart, LookStop};
end