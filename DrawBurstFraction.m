function DrawBurstFraction(Organized)
Conditions = Organized.Conditions;
for n = 1:length(Conditions)
  Condition = Conditions{n};
  if(strcmp(Condition, 'ptx'))
    continue;
  end
  MakeBurstFractionFigure(Organized, Condition);
  MakeBurstFractionFigure(Organized, 'ptx', Condition);
end
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MakeBurstFractionFigure(Organized, Condition, NonControl)
g_syn = 10:15:100;
g_h = 10:15:100;

%find mean control frequency
Temp = zeros(length(g_syn), length(g_h));
Averages = Temp;
TypeArr = Organized.(sprintf('Type_%s', Condition));
if(nargin == 3)
  NonControlArray = Organized.(sprintf('Type_%s', NonControl));
  AdmitArray = isfinite(NonControlArray);
  Title = sprintf('Fraction of Bursters for %s control', NonControl);
else
  AdmitArray = ones(length(Organized.g_h), 1);
  Title = sprintf('Fraction of Bursters for %s', Condition);
end
for n = 1:length(Organized.g_h)
  gs = find(g_syn == Organized.g_syn(n));
  if(length(gs) == 0)
    continue;
  end
  gh = find(g_h == Organized.g_h(n));
  if(length(gh) == 0)
    continue;
  end
  if(~isfinite(TypeArr(n)) | ~AdmitArray(n))
    continue;
  end
  
  Temp(gs, gh) = Temp(gs, gh) + 1;
  if(TypeArr(n) > 0)
    Averages(gs, gh) = Averages(gs, gh) + 1;
  end
end
gInd = find(isfinite(Averages) & Temp > 0 & Averages > 0);
Temp = Averages ./ Temp;
[x, y] = ind2sub(size(Temp), gInd);
x = g_syn(x);
y = g_h(y);
Colors = Temp(gInd);

h = NamedFigure(Title);
set(h, 'WindowStyle', 'docked');
clf;
hold on;
CircleSize = 1000;
scatter(x, y, CircleSize, Colors, 'filled');
Map = repmat((200:-1:0)'/200, 1, 3);
colormap(Map);
colorbar;
caxis([0 1]);

xlim([0 140])
ylim([0 140])
axis square;
xlabel('g_s_y_n');
ylabel('g_h');
title(Title);
Axes_h = get(h, 'CurrentAxes');
set(Axes_h, 'FontName', 'Arial');
set(Axes_h, 'FontSize', [30]);
set(Axes_h, 'YTick', 20:20:140);
set(Axes_h, 'XTick', 0:20:140);
return
