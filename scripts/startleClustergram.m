function [cluster,exp] = startleClustergram(exp,varargin)
p = inputParser;
p.addRequired('exp',@(x) isstruct(x) || exist(x,'dir'));
p.addParameter('lights',{'on','off'},@iscellstr);
p.addParameter('normalize',1,@(x) islogical(x) || isnumeric(x));
p.addParameter('normSearch',{'wt' '+'},@iscellstr);
p.addParameter('overwrite',1,@(x) islogical(x) || isnumeric(x));
p.addParameter('plot',1,@(x) islogical(x) || isnumeric(x));
p.addParameter('includeMutant',{},@iscellstr);
p.addParameter('includeGeno',{},@iscellstr);
p.addParameter('param',{...
	'STIM_mean',...
	'POST_mean',...
	'POST_std',...
	'S1_mean',...
	'P1_mean',...
	'S2_mean',...
	'P2_mean',...
	'S3_mean',...
	'P3_mean',...
	'S4_mean',...
	'P4_mean',...
	'S5_mean',...
	'P5_mean',...
	'STIMvPOST_mean',...
	'BASELINEvsSTIM_mean',...
     },@iscellstr);
p.parse(exp,varargin{:});
inputs = p.Results;

if ~isstruct(exp) && exist(exp,'dir')
    fprintf('Loading preprocessed .mats...\n');
    exp = scrapeExp(exp,'match',true);
    fprintf('\n');
end

% include - prune mutants
if ~isempty(inputs.includeMutant)
    exp = exp(cellfun(@(x) any(strcmp(inputs.includeMutant,x)),{exp.mutant}));
end

% load mat files with animalBehaviorAnalysis objects
cluster.norm.on = {};
cluster.norm.off = {};
cluster.rows = {};
cluster.cols = {};
cluster.data = [];

% rearrange data into exp.geno.cluster struct
for i=1:numel(exp)
    fprintf('%s\t%s',exp(i).mutant,exp(i).name);
    % load .mat file
    tmp = load(exp(i).mat,'n');
    obj = tmp.n;
    % extract stats per genotype (exp.geno.cluster.raw)
    t = obj.statistics.uni;
    if ~isfield(exp(i),'geno') || isempty(exp(i).geno)
        exp(i).geno = obj.getGenotypes;
    end
    for j=1:numel(exp(i).geno)
        row_idx = ismember(t.group_label2,exp(i).geno(j).label);
        col_idx = ismember(t.Properties.VariableNames,inputs.param);
        tmp = t(row_idx,col_idx);
        exp(i).geno(j).cluster.columns = t.Properties.VariableNames(col_idx);
        exp(i).geno(j).cluster.raw = table2array(tmp);
        exp(i).geno(j).cluster.norm = [];
        exp(i).geno(j).cluster.data = [];
    end
    exp(i).geno = fixGeno(exp(i).geno);
    % normalize
    if inputs.normalize
        % find normalization index
        geno = {exp(i).geno.label};
        n = zeros(size(geno));
        for j=1:numel(geno)
            n(j) = max(cellfun(@numel,regexpi(geno{j},inputs.normSearch)));
        end
        [max_n,norm_idx] = max(n);
        if ~any(n)
            fprintf('\tNORM FAILED: no genotype found');
        elseif sum(max_n == n) > 1
            fprintf('\tNORM FAILED: %d genotypes found',sum(max_n==n));
        else
            % get norm vals
            tmp = exp(i).geno(norm_idx).cluster.raw;
            norm.mean = mean(tmp,1);
            norm.std = std(tmp,0,1);
            norm.n = size(tmp,1);
            % normalize each geno
            for j=1:numel(exp(i).geno)
                exp(i).geno(j).cluster.norm = (exp(i).geno(j).cluster.raw - norm.mean ) ./ norm.std;
                exp(i).geno(j).cluster.data = mean(exp(i).geno(j).cluster.norm,1);
                if any(isnan(exp(i).geno(j).cluster.data)) || any(abs(exp(i).geno(j).cluster.data)==inf) 
                    fprintf('\tNORM FAILED: inf/nan values resulting');
                    exp(i).geno(j).cluster.norm = [];
                end
            end
        end
    else
        for j=1:numel(exp(i).geno)
            exp(i).geno(j).cluster.data = mean(exp(i).geno(j).cluster.raw,1);
        end
    end
    fprintf('\n');
end
% setup data
for i=1:numel(exp)
    for j=1:numel(exp(i).geno)
        g = exp(i).geno(j).label;
        l = exp(i).lights;
        % find other indices
        switch l
            case 'on'
                lo = 'off';
            case 'off'
                lo = 'on';
        end
        idx_exp = ismember({exp.name},regexprep(exp(i).name,l,lo)) & ismember({exp.mutant},exp(i).mutant);
        idx_geno = ismember({exp(idx_exp).geno.label},exp(i).geno(j).label);
        idx = find(ismember(cluster.rows,g));
        % add to global output
        if ~isempty(exp(i).geno(j).cluster.norm) && sum(idx_geno) == 1 && ~isempty(exp(idx_exp).geno(idx_geno).cluster.norm)
            if ~isempty(idx)
                if numel(cluster.norm.(l)) < idx
                    cluster.norm.(l){idx} = exp(i).geno(j).cluster.norm;
                else
                    cluster.norm.(l){idx} = [cluster.norm.(l){idx}; exp(i).geno(j).cluster.norm];
                end
            else
                cluster.norm.(l){end+1} = exp(i).geno(j).cluster.norm;
                cluster.rows{end+1} = g;
            end
        end
    end
end
% aggregate data
on = cellfun(@mean,cluster.norm.on,'UniformOutput',0);
off = cellfun(@mean,cluster.norm.off,'UniformOutput',0);
if numel(on) ~= numel(off)
   error('Different number of genotypes in ''on'' and ''off'''); 
elseif range(cellfun(@numel,[on off]))
   error('Diffeent number of parameters in ''on'' and/or ''off'' data'); 
end
for i=1:numel(on)
    cluster.data(i,1:numel(on{1})) = on{i};
    cluster.data(i,numel(on{1})+1:numel(on{1})*2) = off{i};
end
cluster.cols = [cellfun(@(x) [x '_on'],inputs.param,'UniformOutput',0) cellfun(@(x) [x '_off'],inputs.param,'UniformOutput',0)];
% plot
if inputs.plot
    % include
    if ~isempty(inputs.includeGeno)
       idx = ismember(cluster.rows,inputs.includeGeno);
       cluster.rows = cluster.rows(idx);
       cluster.data = cluster.data(idx,:);
       fprintf('%d/%d ''include'' genotypes\n',sum(idx),numel(inputs.includeGeno));
    end
    rm = find(isnan(mean(cluster.data,1)) | abs(mean(cluster.data,1)) == inf);
    for i=1:numel(rm)
       fprintf('%s removed for inf/nan values\n',cluster.cols{rm(i)}); 
    end
    cluster.data(:,rm) = [];
    cluster.cols(rm) = [];
    
    HeatMap(cluster.data,'columnlabels',regexprep(cluster.cols,{'_'},' '),'columnlabelsrotate',45,'rowlabels',cluster.rows,'displayrange',3)
    clustergram(cluster.data,'columnlabels',regexprep(cluster.cols,{'_'},' '),'columnlabelsrotate',45,'rowlabels',cluster.rows,'displayrange',3);
end

end

function exp = scrapeExp(exp_path,varargin)
p = inputParser;
p.addRequired('exp_path',@(x) exist(x,'dir'));
p.addParameter('match',false,@(x) islogical(x) || isnumeric(x));
p.addParameter('lights',{'on','off'},@iscellstr);
p.parse(exp_path,varargin{:});
inputs = p.Results;

% get all mutant dirs
d = dir(exp_path);
d(~[d.isdir]) = [];
d(cellfun(@(x) strcmp(x(1),'.'),{d.name})) = [];

% create exp struct
exp = [];
for i=1:numel(d)
    % get experiments
    d1 = dir(fullfile(d(i).folder,d(i).name));
    d1(~[d1.isdir]) = [];
    d1(cellfun(@(x) strcmp(x(1),'.'),{d1.name})) = [];
    d1(contains({d1.name},'combined')) = [];
    for j=1:numel(d1)
        mat = dir(fullfile(d1(j).folder,d1(j).name,'*.mat'));
        if numel(mat) == 1
            exp(end+1).name = d1(j).name;
            exp(end).mutant = d(i).name;
            exp(end).path = fullfile(d1(j).folder,d1(j).name);
            exp(end).mat = fullfile(exp(end).path,mat.name);
            lights = strsplit(regexprep(mat.name,{'_geno','\.mat'},''),'_');
            lights = lights{end};
            if ~any(strcmp(lights,inputs.lights))
                if strcmp(lights,'err')
                    fprintf('SKIPPED: error during preprocessing');
                else
                    fprintf('ERROR: %s isn''t a light parameter',lights);
                end
                fprintf('\t%s\t%s\n',exp(end).mutant,exp(end).name);
                
                exp(end) = [];
            else
                exp(end).lights = lights;
            end
        else
            fprintf('ERROR: %d .mats\t%s\t%s\n',numel(mat),d(i).name,d1(j).name);
        end
    end
end
% remove non paired on/off
if inputs.match
    rm = false(size(exp));
    for i=1:numel(exp)
       name = regexprep(exp(i).name,['_' exp(i).lights],'');
       switch exp(i).lights
           case 'on'
               lights = 'off';
           case 'off'
               lights = 'on';
       end
       name = [name '_' lights];
       if ~ismember(name,{exp.name})
           rm(i) = true;
           fprintf('REMOVED: %s doesn''t have a lights ''%s'' match\n',exp(i).name,lights);
       end
    end
    exp(rm) = [];
end

end

function geno = fixGeno(geno,varargin)
% inputs
p = inputParser;
p.addParameter('geno',{...
    'HET',...
    'HOM',...
    'WT',...
	'HET:HET',...
	'HET:HOM',...
	'HET:WT',...
	'HOM:HET',...
	'HOM:WT',...
	'HOM:HOM',...
	'WT:HET',...
	'WT:HOM',...
	'WT:WT'},...
    @iscellstr);
p.parse(varargin{:});
inputs = p.Results;

% get labels and fix common typos
label = {geno.label};
label = cellfun(@strtrim,label,'UniformOutput',0);
label = regexprep(label,{': ',' :'},':');
label = regexprep(label,{'/ ',' /'},'/');
label = regexprep(label,'del ','del');
label = regexprep(label,'SCN','scn ');
label = regexprep(label,'Scn1lab','scn');
label = regexprep(label,'scn del44/\+','scn del44 HET');
label = regexprep(label,'scn del5/\+','scn del5 HET');

% check if can be fixed (all genotypes in isolation with one preceded by gene)
idx = ismember(label,inputs.geno);
if sum(idx) == numel(idx) - 1
   % split prefix from genotype
   idx = find(~idx);
   tmp = strsplit(label{idx});
   prefix = lower(join(tmp(1:end-1)));
   prefix = prefix{1};
   label{idx} = tmp{end};
   for i=1:numel(label)
      label{i} = [prefix ' ' label{i}];
   end
end
[geno.label] = label{:};
end


