function PlotPooled(Organized)
tic
ValueString = 'HalfCenter.Freq';  %Tells which values to map
if(nargin < 1)
  OrganizedExists = evalin('base', 'exist(''Organized'')');
  if(OrganizedExists)
    Organized = evalin('base', 'Organized');
  else
    Organized = OrganizeByParameters([], ValueString);
    assignin('base', 'Organized', Organized);
  end
end

close all

DrawMeanChange(Organized);
DrawScatter(Organized);
DrawBurstFraction(Organized);
DrawBar(Organized);

Props = LoadIntrinsicProperties;
DrawCorrelateProps(Organized, Props);
return
toc