% script to create animalBehaviorAnalysis object and perform built-in analyses 
%
% SurveyBott, 2018, info@surveybott.com

function cmd = run_aba(cmd,varargin)
% inputs
p = inputParser;
p.addRequired('cmd',@(x) isstruct(x) && all(isfield(x,{'mutant','lights','files','savePath','projectName'})) && all(isfield(x.mutant,{'name'})));
p.addParameter('acclimationTiming',[0 4800],@(x) isnumeric(x) && numel(x) == 2);
p.addParameter('baselineTiming',[4800 5400],@(x) isnumeric(x) && numel(x) == 2);
p.addParameter('eventTiming',[0 1 29],@(x) isnumeric(x) && numel(x) == 3);
p.addParameter('eventIndex','Light',@ischar);
p.addParameter('dv','actinteg',@ischar);
p.addParameter('stats',false,@(x) isnumeric(x) || islogical(x));
p.addParameter('moveError',true,@(x) isnumeric(x) || islogical(x));
p.addParameter('normalize',true,@(x) isnumeric(x) || islogical(x));
p.addParameter('liteTs',true,@(x) isnumeric(x) || islogical(x));
p.parse(cmd,varargin{:});
inputs = p.Results;

% save location
projectName = cmd.projectName;
saveLocation = fullfile(cmd.savePath);
if ~exist(saveLocation,'dir')
    % handle permission issuess
    try
        mkdir(saveLocation);
    catch
        saveLocation = [saveLocation '_new'];
        mkdir(saveLocation);
        cmd.savePath = saveLocation;
    end
end

% setup aba inputs
cmd.err = [];
try
if isstruct(cmd.files)
    % load
    groupingFile = cmd.files.grouping;
    viewPointFile = cmd.files.viewpoint;
    timingFile = cmd.files.timing;
    if iscell(timingFile)
        timingFile = timingFile{1};
    end
    eventLabel = cmd.lights;
    n = animalBehaviorAnalysis('PROJECT_NAME',projectName,...
        'SAVE_LOCATION',saveLocation,...
        'ACCLIMATION_TIMING', inputs.acclimationTiming,...
        'EVENT_TIMING',inputs.eventTiming,...
        'BASELINE_TIMING',inputs.baselineTiming,...
        'GROUPING_FILE', groupingFile,...
        'VIEWPOINT_FILE',viewPointFile,...
        'TIMING_FILE',timingFile,...
        'EVENT_INDEX',inputs.eventIndex,...
        'EVENT_LABEL',eventLabel,...
        'DEP_VAR',inputs.dv);
elseif iscellstr(cmd.files)
    % merge .mat
    n = animalBehaviorAnalysis('PROJECT_NAME',projectName,'SAVE_LOCATION',saveLocation,'MERGE',cmd.files);
end
if inputs.stats
    n.computeAnimalWiseProperties;
    if inputs.liteTs
        fig = fullfile(saveLocation,projectName,'results','Full_Raw_Timeseries.fig');
        if exist(fig,'file')
            liteTimeSeries(fig);
        end
        fig = fullfile(saveLocation,projectName,'results','NoBaseline_Raw_Timeseries.fig');
        if exist(fig,'file')
            liteTimeSeries(fig);
        end
    end
end

catch err
    cmd.err = err;
end

% strip down animalBehaviorAnalysis object
n.data.viewpoint_data = [];

% save .mat of animalBehaviorAnalysis obj and cmd struct
saveFile = fullfile(saveLocation,projectName,projectName);
if ~isempty(cmd.err)
    saveFile = [saveFile '_err'];
    if inputs.moveError && isstruct(cmd.files) 
       dataPath = cmd.mutant.path;
       [errPath,mutantFolder] = fileparts(dataPath);
       errPath = fullfile(fileparts(errPath),'error',mutantFolder);
       dataPath = cmd.exp.path;
       if ~isdir(errPath)
           mkdir(errPath);
           system(sprintf('chmod 770 "%s"',errPath));
       end
       % copy
       code = system(sprintf('rsync -av "%s" "%s"',dataPath,errPath));
       % open permissions
       system(sprintf('chmod -R 770 "%s"',fullfile(errPath)));
       % delete
       if ~code
           system(sprintf('rm -rf "%s"',dataPath));
           try
               rmdir(fileparts(dataPath));
           catch
           end
       end
    end
elseif numel(cmd.exp) == 1
   cmd.exp.analysis = [saveFile '.mat']; 
end

saveFile = [saveFile '.mat'];
try
    save(saveFile,'n','cmd');
catch
    newSaveFile = regexprep(saveFile,'\.mat','');
    if ~strcmp(newSaveFile(end-2:end),'err')
        newSaveFile = [newSaveFile '_err.mat'];
        try
            movefile(saveFile,newSaveFile);
        catch
        end
    end
end

% remove err.mat, if it exists
if ~contains(saveFile,'_err.mat')
   errFile = regexprep(saveFile,'\.mat','_err\.mat');
   try
      delete(errFile); 
   catch
   end
end

end