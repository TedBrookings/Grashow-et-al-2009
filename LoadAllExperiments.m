function Experiments = LoadAllExperiments(ListFile, BaseDir)
%Experiments = LoadAllExperiments(ListFile, BaseDir)
% returns a list of experiment data
%    -BaseDir is the directory for all the experiment foloders
%         (defaults to C:\Documents and Settings\Rachel\My Documents\spectra)
%    -ListFile is a text file listing all the directories to open
%         (defaults to folder_names.txt)

if(nargin > 3)
  error('LoadAllExperiments.m requires 2 or fewer input arguments');
end
if(nargin < 2)
  if(ispc)
    BaseDir = 'C:\Documents and Settings\Rachel\My Documents\spectra\';
    AltDir = '/mnt/dwidget/';
  else
    BaseDir = '/mnt/dwidget/';
    AltDir = 'C:\Documents and Settings\Rachel\My Documents\spectra\';
  end
end
if(nargin < 1)
  ListFile = [BaseDir, 'folder_names.txt'];
end


[DirList, ExpList] = GetDirList(BaseDir, ListFile);
Experiments = GetExperimentList(DirList, ExpList, BaseDir, AltDir);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [DirList, ExpList] = GetDirList(BaseDir, ListFile)
DirList = {};
ExpList = {};
fid = fopen(ListFile, 'r');
NextLine = fgetl(fid);
while ischar(NextLine);
  ExpList = {ExpList{:}, NextLine};
  DirList = {DirList{:}, [BaseDir, NextLine]};
  NextLine = fgetl(fid);
end
fclose(fid);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Experiments = GetExperimentList(DirList, ExpList, BaseDir, ...
					 AltDir)
Experiments = [];
if(ispc)
  Slash = '\';
  AltSlash = '/';
else
  Slash = '/';
  AltSlash = '\';
end

for n = 1:length(DirList)
  FileName = [DirList{n}, Slash, 'Analysis.mat'];
  try
    load(FileName);
    for m=1:length(Analysis)
      FindInd = strfind(Analysis(m).FileName, AltSlash);
      if(length(FindInd) > 1)
	FindInd = FindInd(end);
	NewFile = Analysis(m).FileName;
	NewFile = NewFile((FindInd+1):end);
	NewFile = [BaseDir, NewFile];
	Analysis(m).FileName = NewFile;
      end
    end
    Experiment.Analysis = Analysis;
    Experiment.ExpNum = ExpNum;
    Experiments = [Experiments, Experiment];
  catch
    WarnStr = sprintf('Warning, %s does not exist, so data from %s will not be used!', ...
                        FileName, ExpList{n});
    disp(WarnStr)
  end
end
return