function varargout = StencilDeriv(Y, DeltaX)
% [Y', Y'', Y''', Y''''] = StencilDeriv(Y, DeltaX)
%   Uses five point stencil to take derivatives of Y, up to order 4.
%   X is assumed to have constant spacing DeltaX.
% EXPERIMENTAL:
% [Y', Y''] = StencilDeriv(Y, X)
%   Uses Lagrange Interpolating Polynomial to take derivatives up
%   to order 2, when X is not assumed to be evenly spaced.  Not
%   very stable to variation in spacing of X.
if(nargout < 1 | nargout > 4)
  error('Incorrect number of output arguments.  Run "help StencilDeriv".')
end
if(nargin ~= 2)
  error('Incorrect number of input arguments.  Run "help StencilDeriv".')
end

if(length(DeltaX) > 1)
  varargout = VariableStencil(Y, DeltaX, nargout);
  return
end
  
if(size(Y,2) > size(Y,1))
  Y = Y';
  Flip = true;
else
  Flip = false;
end

if(nargout == 1)
  YPrime = zeros(size(Y));
  YPrime(3:(end-2)) = ( Y(1:(end-4),:) - Y(5:end,:) + ...
		8 * (Y(4:(end-1),:) - Y(2:(end-3),:)) ) / (12 * DeltaX);
  YPrime(1) = (-25*Y(1) + 48*Y(2) - 36*Y(3) + 16*Y(4) - 3*Y(5)) ...
      / (12 * DeltaX);
  YPrime(2) = (-3*Y(1) - 10*Y(2) + 18*Y(3) - 6*Y(4) + Y(5)) ...
      / (12 * DeltaX);
  YPrime(end-1) = (3*Y(end) + 10*Y(end-1) - 18*Y(end-2) + 6*Y(end-3) ...
		   - Y(end-4)) / (12 * DeltaX);
  YPrime(end) = (25*Y(end) - 48*Y(end-1) + 36*Y(end-2) - 16*Y(end-3) ...
		 + 3*Y(end-4)) / (12 * DeltaX);
  
  varargout = {YPrime};
else
  Y_Hm2 = Y(1:(end-4),:);
  Y_Hm1 = Y(2:(end-3),:);
  Y_H = Y(3:(end-2),:);
  Y_Hp1 = Y(4:(end-1),:);
  Y_Hp2 = Y(5:end,:);
  
  YSize = size(Y);
  Y_d1 = zeros(YSize);
  Y_d1(3:(end-2),:) = (Y_Hm2 - Y_Hp2 + 8 * (Y_Hp1 - Y_Hm1)) / (12 * DeltaX);
  Y_d1(1,:) = (-25*Y(1) + 48*Y(2) - 36*Y(3) + 16*Y(4) - 3*Y(5)) ...
	    / (12 * DeltaX);
  Y_d1(2,:) = (-3*Y(1) - 10*Y(2) + 18*Y(3) - 6*Y(4) + Y(5)) ...
	    / (12 * DeltaX);
  Y_d1(end-1,:) = (3*Y(end) + 10*Y(end-1) - 18*Y(end-2) + 6*Y(end-3) ...
		 - Y(end-4)) / (12 * DeltaX);
  Y_d1(end,:) = (25*Y(end) - 48*Y(end-1) + 36*Y(end-2) - 16*Y(end-3) ...
	       + 3*Y(end-4)) / (12 * DeltaX);
  
  DX2 = DeltaX * DeltaX;
  Y_d2 = zeros(YSize);
  Y_d2(3:(end-2),:) = (16 * (Y_Hm1 + Y_Hp1) - (Y_Hm2 + Y_Hp2) - 30 * Y_H) ...
      / (12 * DX2);
  Y_d2(1,:) = (35*Y(1) - 104*Y(2) + 114*Y(3) - 56*Y(4) + 11*Y(5)) ...
	    / (12 * DX2);
  Y_d2(2,:) = (11*Y(1) - 20*Y(2) + 6*Y(3) + 4*Y(4) - Y(5)) ...
	    / (12 * DX2);
  Y_d2(end-1,:) = (11*Y(end) - 20*Y(end-1) + 6*Y(end-2) + 4*Y(end-3) ...
		 - Y(end-4)) / (12 * DX2);
  Y_d2(end,:) = (35*Y(end) - 104*Y(end-1) + 114*Y(end-2) - 56*Y(end-3) ...
	       + 11*Y(end-4)) / (12 * DX2);
  if(nargout == 2)
    varargout = {Y_d1, Y_d2};
  else
    DX3 = DX2 * DeltaX;
    Y_d3 = zeros(YSize);
    Y_d3(3:(end-2),:) = (Y_Hp2 - Y_Hm2 + 2 * (Y_Hm1 - Y_Hp1)) / (2 * DX3);
    Y_d2(1,:) = (-5*Y(1) + 18*Y(2) - 24*Y(3) + 14*Y(4) - 3*Y(5)) ...
	      / (2 * DX3);
    Y_d2(2,:) = (-3*Y(1) + 10*Y(2) - 12*Y(3) + 6*Y(4) - Y(5)) ...
	      / (2 * DX3);
    Y_d2(end-1,:) = (3*Y(end) - 10*Y(end-1) + 12*Y(end-2) - 6*Y(end-3) ...
		   + Y(end-4)) / (2 * DX3);
    Y_d2(end,:) = (5*Y(end) - 18*Y(end-1) + 24*Y(end-2) - 14*Y(end-3) ...
		 + 3*Y(end-4)) / (2 * DX3);
    if(nargout == 3)
      varargout = {Y_d1, Y_d2, Y_d3};
    else
      DX4 = DeltaX * DX3;
      Y_d4 = zeros(YSize);
      Y_d4(3:(end-2),:) = (Y_Hp2 + Y_Hm2 - 4 * (Y_Hm1 + Y_Hp1) + 6 * Y_H)/DX4;
      Y_d4(1,:) = Y_d4(3);
      Y_d4(2,:) = Y_d4(3);
      Y_d4(end-1,:) = Y_d4(end-2);
      Y_d4(end,:) = Y_d4(end-1);
      varargout = {Y_d1, Y_d2, Y_d3, Y_d4};
    end
  end
end
if(Flip)
  for n=1:nargout
    varargout{n} = varargout{n}';
  end
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Output = VariableStencil(Y, X, NumArgOut)

NumX = length(X);
if(NumX ~= length(Y))
  error('X and Y must have the same number of sample points')
end

Y_d1 = zeros(size(Y));
if(NumArgOut > 1)
  Y_d2 = zeros(size(Y));
end
X1 = X(1); X2 = X(2); X3 = X(3); X4 = X(4); X5 = X(5);
Y1 = Y(1); Y2 = Y(2); Y3 = Y(3); Y4 = Y(4); Y5 = Y(5);
D1 = (X2-X1)*(X3-X1)*(X4-X1)*(X5-X1);
D2 = (X2-X1)*(X3-X2)*(X4-X2)*(X5-X2);
D3 = (X3-X1)*(X3-X2)*(X4-X3)*(X5-X3);
D4 = (X4-X1)*(X4-X2)*(X4-X3)*(X5-X4);
D5 = (X5-X1)*(X5-X2)*(X5-X3)*(X5-X4);

%n=1:
C1 = -1.0/(X2-X1) - 1.0/(X3-X1) - 1.0/(X4-X1) - 1.0/(X5-X1);
C2 = (X3-X1)*(X4-X1)*(X5-X1) / D2;
C3 = -(X2-X1)*(X4-X1)*(X5-X1) / D3;
C4 = (X2-X1)*(X3-X1)*(X5-X1) / D4;
C5 = -(X2-X1)*(X3-X1)*(X4-X1) / D5;
Y_d1(1) = C1*Y1 + C2*Y2 + C3*Y3 + C4*Y4 + C5*Y5;


%n=2:
C1 = -(X3-X2)*(X4-X2)*(X5-X2)/ D1;
C2 = 1.0/(X2-X1) - 1.0/(X3-X2) - 1.0/(X4-X2) - 1.0/(X5-X2);
C3 = (X2-X1)*(X4-X2)*(X5-X2)/ D3;
C4 = -(X2-X1)*(X3-X2)*(X5-X2)/ D4;
C5 = (X2-X1)*(X3-X2)*(X4-X2)/ D5;
Y_d1(2) = C1*Y1 + C2*Y2 + C3*Y3 + C4*Y4 + C5*Y5;

for n = 3:(NumX-2)
  if(n > 3)
    X1 = X2; X2 = X3; X3 = X4; X4 = X5; X5 = X(n + 2);
    Y1 = Y2; Y2 = Y3; Y3 = Y4; Y4 = Y5; Y5 = Y(n + 2);
    D1 = (X2-X1)*(X3-X1)*(X4-X1)*(X5-X1);
    D2 = (X2-X1)*(X3-X2)*(X4-X2)*(X5-X2);
    D3 = (X3-X1)*(X3-X2)*(X4-X3)*(X5-X3);
    D4 = (X4-X1)*(X4-X2)*(X4-X3)*(X5-X4);
    D5 = (X5-X1)*(X5-X2)*(X5-X3)*(X5-X4);
  end
  
  C1 = (X3-X2)*(X4-X3)*(X5-X3) / D1;
  C2 = -(X3-X1)*(X4-X3)*(X5-X3) / D2;
  C3 = 1.0/(X3-X1) + 1.0/(X3-X2) - 1.0/(X4-X3) - 1.0/(X5-X3);
  C4 = (X3-X1)*(X3-X2)*(X5-X3) / D4;
  C5 = -(X3-X1)*(X3-X2)*(X4-X3) / D5;
  Y_d1(n) = C1*Y1 + C2*Y2 + C3*Y3 + C4*Y4 + C5*Y5;
  if(NumArgOut > 1)
    C1 = 2*(3*X3*X3 - 2*X3*(X2+X4+X5) + X2*(X4+X5) + X4*X5) / D1;
    C2 = -2*(3*X3*X3 - 2*X3*(X1+X4+X5) + X1*(X4+X5) + X4*X5) / D2;
    C3 = 2*(-3*X3*(X1+X2+X4+X5-2*X3) + (X1+X2)*(X4+X5) + X1*X2 ...
	    + X4*X5) / D3;
    C4 = -2*(3*X3*X3 - 2*X3*(X1+X2+X5) + X5*(X1+X2) + X1*X2) / D4;
    C5 = 2*(3*X3*X3 - 2*X3*(X1+X2+X4) + X4*(X1+X2) + X1*X2) / D5;
    Y_d2(n) = C1*Y1 + C2*Y2 + C3*Y3 + C4*Y4 + C5*Y5;
  end
end

%n=NumX-1:
C1 = -(X4-X2)*(X4-X3)*(X5-X4) / D1;
C2 = (X4-X1)*(X4-X3)*(X5-X4) / D2;
C3 = -(X4-X1)*(X4-X2)*(X5-X4) / D3;
C4 = 1.0/(X4-X1) + 1.0/(X4-X2) + 1.0/(X4-X3) - 1.0/(X5-X4);
C5 = (X4-X1)*(X4-X2)*(X4-X3) / D5;
Y_d1(NumX - 1) = C1*Y1 + C2*Y2 + C3*Y3 + C4*Y4 + C5*Y5;

%n=NumX:
C1 = (X5-X2)*(X5-X3)*(X5-X4) / D1;
C2 = -(X5-X1)*(X5-X3)*(X5-X4) / D2;
C3 = (X5-X1)*(X5-X2)*(X5-X4) / D3;
C4 = -(X5-X1)*(X5-X2)*(X5-X3) / D4;
C5 = 1.0/(X5-X1) + 1.0/(X5-X2) + 1.0/(X5-X3) + 1.0/(X5-X4);
Y_d1(NumX) = C1*Y1 + C2*Y2 + C3*Y3 + C4*Y4 + C5*Y5;

if(NumArgOut == 1)
  Output = {Y_d1};
elseif(NumArgOut == 2)
  Output = {Y_d1, Y_d2};
end
return