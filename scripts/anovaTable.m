function t = anovaTable(in,varargin)
% inputs
p = inputParser;
p.addRequired('in',@isdir);
p.addParameter('include',{},@iscellstr);
p.parse(in,varargin{:});
inputs = p.Results;

% get .mat files
d = dir(fullfile(in,'**/*.mat'));
d(contains({d.name},'combined')) = [];
d = d(contains({d.name},'_on.mat') | contains({d.name},'_off.mat'));
if ~isempty(inputs.include)
   d = d(contains({d.name},inputs.include) | contains({d.folder},inputs.include));
end

% prepopulate cell array
load(fullfile(d(1).folder,d(1).name));
f = ['folder' 'experiment' 'lights' fieldnames(n.results.oneway)' fieldnames(n.results.twoway)'];
c = cell(numel(d),numel(f));
empty = false(size(d));

% get p values from experiments
parfor i=1:size(c,1)
    % get meta info
    fprintf('%d',i);
    [~,folder] = fileparts(fileparts(d(i).folder));
    [~,tmp] = fileparts(d(i).name);
    r = regexp(tmp,'_');
    experiment = tmp(1:r(end)-1);
    lights = tmp(r(end)+1:end);
    % load file, get p-values
    tmp = num2cell(nan(size(c(i,:))));
    mat = load(fullfile(d(i).folder,d(i).name));
    n = mat.n;
    anova = {'oneway','twoway'};
    for j=1:numel(anova)
        a = anova{j};
        try
            fa = fieldnames(n.results.(a));
            for k=1:numel(fa)
                try
                    switch a
                        case 'oneway'
                            pval = n.results.(a).(fa{k}).p;
                        case 'twoway'
                            pval = n.results.(a).(fa{k}).within_tbl(2,:).pValue;
                    end
                    tmp{ismember(f,fa{k})} = pval; % put into correct place based on header
                catch
                end
            end
        catch
        end
    end
    % add meta info if p values were found
    if ~all(cellfun(@isnan,tmp))
        tmp{1} = folder;
        tmp{2} = experiment;
        tmp{3} = lights;
    else
       empty(i) = true; 
    end
    c(i,:) = tmp;
    fprintf('\n');
end

t = cell2table(c(~empty,:));
t.Properties.VariableNames = f;
end