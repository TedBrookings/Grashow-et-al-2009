function AnalyzeAllExperiments
%  data in folder_names.txt has been moved, so can't analyze
%AnalyzeFromFolderList('folder_names.txt', false);

AnalyzeFromFolderList('MorrisLecarFolders.txt', false);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function AnalyzeFromFolderList(DirListFile, ForceAnalyze)
if(nargin < 2)
  ForceAnalyze = false;
end

if(ispc)
  Slash = '\';
  BaseReadDir = 'C:\Documents and Settings\Rachel\My Documents\spectra';
  BaseWriteDir = 'C:\Documents and Settings\Rachel\My Documents\spectra';
else
  Slash = '/';
  BaseReadDir = '/mnt/dwidget';
  BaseWriteDir = '/mnt/dwidget';
  %  BaseWriteDir = '/home/ted/data/bizarro/';
end

[FolderList, ExpList] = GetFolderList(BaseReadDir, DirListFile);

for n = 1:length(FolderList)
  WriteDir = [BaseWriteDir, Slash, ExpList{n}];
  DirName = FolderList{n};
  [Success, Message, MessageID] = mkdir(WriteDir);
  OutFileName = [WriteDir, Slash, 'Analysis.mat'];
  
  if(~ForceAnalyze & AlreadyDone(OutFileName))
    disp(sprintf('########### Previously Completed %s ###########', DirName))
    continue;
  end
  disp(sprintf('################ Analyzing %s #################', DirName))
  try
    [Analysis, ExpNum] = AnalyzeExperiment(DirName);
  catch
    disp(sprintf('######## Warning... Unable to finish %s ########', DirName))
    ErrStruct = lasterror;
    fprintf(2, [ErrStruct.message, '\n'])
    ErrLoc = ErrStruct.stack(1);
    ErrLine = sprintf('in %s line %g\n', ErrLoc.file, ErrLoc.line);
    fprintf(2, ErrLine)
    continue;
  end

  try
    save(OutFileName, 'Analysis', 'ExpNum');
    disp(sprintf('########### Finished %s ###########', DirName))
  catch
    disp(sprintf('######## Warning... Unable to save %s ########', DirName))
  end
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [FolderList, ExpList] = GetFolderList(BaseDir, DirListFile)
if(ispc)
  Slash = '\';
else
  Slash = '/';
end
DirListFile = [BaseDir, Slash, DirListFile];
FolderList = {};
ExpList = {};
fid = fopen(DirListFile);
TextLine = fgetl(fid);
while(ischar(TextLine))
  FolderList = {FolderList{:}, [BaseDir, Slash, TextLine]};
  ExpList = {ExpList{:}, TextLine};
  TextLine = fgetl(fid);
end
fclose(fid);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Done = AlreadyDone(OutFileName)
fid = fopen(OutFileName, 'r');
Done = (fid > 0);
if(Done)
  fclose(fid);
end
return
