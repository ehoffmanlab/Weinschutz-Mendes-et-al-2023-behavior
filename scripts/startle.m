classdef startle < handle
    properties
        experiment
        day
    end
    methods
        function obj = startle(in,varargin)
            % inputs
            if ~exist('in','var')
                in = pwd;
            end
            p = inputParser;
            p.addRequired('in',@(x) iscellstr(x) || exist(x,'dir'));
            p.addParameter('load',false,@(x) islogical(x) || isnumeric(x));
            p.addParameter('match',false, @(x) islogical(x) || isnumeric(x));
            p.addParameter('includeMutant',{});
            p.parse(in,varargin{:});
            inputs = p.Results;
            % load
            obj.experiment = obj.scrapeExp(in,'match',inputs.match,'includeMutant',inputs.includeMutant);
            if isempty(obj.experiment)
                error('No experiments loaded');
            elseif inputs.load
                obj.load;
            end
        end
        % load mat files with animalBehaviorAnalysis objects, add exp.geno struct with exp.geno.param
        function load(obj,varargin)
            p = inputParser;
            p.addParameter('normalize',true,@(x) islogical(x) || isnumeric(x));
            p.addParameter('normSearch',{'WT' '+'},@iscellstr);
            p.addParameter('normType','z');
            p.addParameter('normOutliers',[]);
            p.addParameter('fixGeno',true,@(x) islogical(x) || isnumeric(x));
            p.addParameter('suffix',true,@(x) islogical(x) || isnumeric(x));
            p.addParameter('rename',{},@iscell);
            p.addParameter('skipLoaded',true,@(x) islogical(x) || isnumeric(x));
            p.addParameter('param',{...
                'STIM_mean',...
                'STIM_std',...
                'STIM_propmv',...
                'STIM_nrmMean',... % BASELINEmSTIM_mean
                'POST_mean',...
                'POST_std',...
                'POST_propmv',...
                'STIM1_mean',...
                'STIM1_propmv',...
                'POST1_mean',...
                'POST1_propmv',...
                'STIM2_mean',...
                'STIM2_propmv',...
                'POST2_mean',...
                'POST2_propmv',...
                'STIM3_mean',...
                'STIM3_propmv',...
                'POST3_mean',...
                'POST3_propmv',...
                'STIM4_mean',...
                'STIM4_propmv',...
                'POST4_mean',...
                'POST4_propmv',...
                'STIM5_mean',...
                'STIM5_propmv',...
                'POST5_mean',...
                'POST5_propmv',...
                'STIMmPOST_mean',...
                'S1mS5_mean',...
                },@iscellstr);
            p.parse(varargin{:});
            inputs = p.Results;
            
            % for each experiment
            for i=1:numel(obj.experiment)
                fprintf('%d/%d\t%s\t%s',i,numel(obj.experiment),obj.experiment(i).mutant,obj.experiment(i).name);
                if inputs.skipLoaded && isfield(obj.experiment,'geno') && isstruct(obj.experiment(i).geno) && isfield(obj.experiment(i).geno(end),'param') && ...
                        isstruct(obj.experiment(i).geno(end)) && isfield(obj.experiment(i).geno(end).param,'data') && ~isempty(obj.experiment(i).geno(end).param.data)
                    fprintf('\talready loaded');
                else
                    % load .mat file, pull indiv data, get groups
                    tmp = load(obj.experiment(i).mat,'n');
                    aba = tmp.n;
                    [~,t] = aba.pullData('SCOPE','INDIV');
                    if ~isfield(obj.experiment(i),'geno') || isempty(obj.experiment(i).geno)
                        label = unique(t.groupLabel);
                        [obj.experiment(i).geno(1:numel(label)).label] = label{:};
                    end
                    % extract parameters per genotype (exp.geno.param.raw)
                    if inputs.fixGeno
                        t.groupLabel = obj.fixGeno(t.groupLabel,'rename',inputs.rename,'suffix',inputs.suffix);
                        label = unique(t.groupLabel);
                        [obj.experiment(i).geno(1:numel(label)).label] = label{:};
                    end
                    for j=1:numel(obj.experiment(i).geno)
                        obj.experiment(i).geno(j).param.columns = {};
                        obj.experiment(i).geno(j).param.raw = [];
                        obj.experiment(i).geno(j).param.norm = [];
                        obj.experiment(i).geno(j).param.data = [];
                        row_idx = ismember(t.groupLabel,obj.experiment(i).geno(j).label);
                        col_idx = ismember(t.Properties.VariableNames,inputs.param);
                        obj.experiment(i).geno(j).animal_id = t.animal_id(row_idx);
                        tmp = t(row_idx,col_idx);
                        if ~isempty(tmp)
                            obj.experiment(i).geno(j).param.columns = t.Properties.VariableNames(col_idx);
                            obj.experiment(i).geno(j).param.raw = table2array(tmp);
                        else
                            fprintf('\n\t\tNO DATA: %s',obj.experiment(i).geno(j).label);
                        end
                    end
                    % normalize
                    if inputs.normalize
                        obj.experiment(i).geno = obj.normGeno(obj.experiment(i).geno,'normSearch',inputs.normSearch,'normType',inputs.normType,'normOutliers',inputs.normOutliers);
                    else
                        for j=1:numel(obj.experiment(i).geno)
                            obj.experiment(i).geno(j).param.data = nanmean(obj.experiment(i).geno(j).param.raw,1);
                            obj.experiment(i).geno(j).param.norm = obj.experiment(i).geno(j).param.data;
                        end
                    end
                end
                fprintf('\n');
            end
        end
        function normalize(obj,varargin)
           p = inputParser;
           p.addParameter('normSearch',{'WT' '+'},@iscellstr);
           p.addParameter('normType','z');
           p.addParameter('normOutliers',[]);
           p.parse(varargin{:});
           inputs = p.Results;
           if obj.isloaded
              for i=1:numel(obj.experiment) 
                 fprintf('(%d/%d)\t%s\t%s',i,numel(obj.experiment),obj.experiment(i).mutant,obj.experiment(i).name);
                 obj.experiment(i).geno = obj.normGeno(obj.experiment(i).geno,...
                     'normSearch',inputs.normSearch,'normOutliers',inputs.normOutliers);
                 fprintf('\n');
              end
           end
        end
        function combineDay(obj,varargin)
            p = inputParser;
            p.addParameter('normalize',true,@(x) islogical(x) || isnumeric(x));
            p.addParameter('normSearch',{'WT' '+'},@iscellstr);
            p.addParameter('normType','z');
            p.addParameter('normOutliers',3);
            p.addParameter('normUpdate',true);
            p.parse(varargin{:});
            inputs = p.Results;
            
            if obj.isloaded
                obj.day = [];
                % get/loop over unique lights
                lights = unique({obj.experiment.lights});
                for x=1:numel(lights)
                    light = obj.experiment(strcmp({obj.experiment.lights},lights{x}));
                    % get/loop over unique mutants
                    m = unique({light.mutant});
                    for i=1:numel(m)
                        % get/loop over unique days
                        d = unique({light(strcmp({light.mutant},m{i})).day});
                        for j=1:numel(d)
                            fprintf('%s\t%s_%s',m{i},d{j},lights{x});
                            d_idx = numel(obj.day) + 1;
                            obj.day(d_idx).mutant = m{i};
                            obj.day(d_idx).day = d{j};
                            obj.day(d_idx).plate = {};
                            obj.day(d_idx).mat = {};
                            obj.day(d_idx).lights = lights{x};
                            obj.day(d_idx).geno = [];
                            % get/loop over experiments in day
                            day = light(strcmp({light.mutant},m{i}) & strcmp({light.day},d{j}));
                            for k=1:numel(day)
                                if k==1
                                    fprintf('\t');
                                end
                                tmp = strsplit(day(k).plate,'_');
                                fprintf('%s',tmp{end});
                                if k < numel(day)
                                    fprintf(',');
                                end
                            end
                            for k=1:numel(day)
                                obj.day(d_idx).plate{end+1} = day(k).plate;
                                obj.day(d_idx).mat{end+1} = day(k).mat;
                                % loop over genotype labels
                                for l=1:numel(day(k).geno)
                                    % find geno idx, add data
                                    if isempty(obj.day(d_idx).geno) || ~sum(strcmp({obj.day(d_idx).geno.label},day(k).geno(l).label))
                                        % new geno
                                        g_idx = numel(obj.day(d_idx).geno) + 1;
                                        geno = day(k).geno(l);
                                        geno.plate = repmat({day(k).plate},numel(geno.animal_id),1);
                                        try
                                            geno.param = rmfield(geno.param,{'data','norm'});
                                        end
                                        if g_idx == 1
                                            obj.day(d_idx).geno = geno;
                                        else
                                            obj.day(d_idx).geno(g_idx) = geno;
                                        end
                                    else
                                        % existing geno
                                        g_idx = strcmp({obj.day(d_idx).geno.label},day(k).geno(l).label);
                                        obj.day(d_idx).geno(g_idx).animal_id = [obj.day(d_idx).geno(g_idx).animal_id; day(k).geno(l).animal_id];
                                        obj.day(d_idx).geno(g_idx).plate = [obj.day(d_idx).geno(g_idx).plate; repmat({day(k).plate},numel(day(k).geno(l).animal_id),1)];
                                        obj.day(d_idx).geno(g_idx).param.raw = [obj.day(d_idx).geno(g_idx).param.raw; day(k).geno(l).param.raw];
                                    end
                                end
                            end
                            % normalize
                            if inputs.normalize
                                obj.day(d_idx).geno = obj.normGeno(obj.day(d_idx).geno,'normSearch',inputs.normSearch,'normType',inputs.normType,...
                                    'normOutliers',inputs.normOutliers,'normUpdate',inputs.normUpdate);
                            else
                                for k=1:numel(obj.day(d_idx).geno)
                                    obj.day(d_idx).geno(k).param.data = nanmean(obj.day(d_idx).geno(k).param.raw,1);
                                    obj.day(d_idx).geno(k).param.norm = obj.day(d_idx).geno(k).param.data;
                                end
                            end
                            fprintf('\n');
                        end
                    end
                end
            else
                error('Experiment data not loaded');
            end
        end
        function bool = isloaded(obj)
            if ~isempty(obj.experiment)  && isfield(obj.experiment,'geno')
                bool = true;
            else
                bool = false;
            end
        end
        function rename(obj,rename,varargin)
            p = inputParser;
            p.addRequired('rename',@iscell);
            p.addParameter('print',true,@(x) islogical(x) || isnumeric(x));
            p.parse(rename,varargin{:});
            inputs = p.Results;
            if obj.isloaded
                for i=1:numel(obj.experiment)
                    if isstruct(obj.experiment(i).geno)
                        obj.experiment(i).geno = obj.renameGeno(obj.experiment(i).geno,rename);
                    end
                end
                if inputs.print
                    obj.printExperiments;
                end
            end
        end
        function printExperiments(obj,varargin)
            p = inputParser;
            p.addParameter('n',true, @(x) islogical(x) || isnumeric(x));
            p.parse(varargin{:});
            inputs = p.Results;
            if ~isempty(obj.experiment)
                e = obj.experiment;
                for i=1:numel(e)
                    fprintf('%s\t%s\n',e(i).mutant,e(i).name)
                    if isfield(e(i),'geno') && isstruct(e(i).geno) && isfield(e(i).geno,'label')
                        for j=1:numel(e(i).geno)
                            fprintf('\t%s',e(i).geno(j).label);
                            if inputs.n && isfield(e(i).geno,'animal_id') && isnumeric(e(i).geno(j).animal_id)
                                fprintf('\tn=%d',numel(e(i).geno(j).animal_id));
                            end
                            fprintf('\n');
                        end
                    end
                end
            end
        end
        function geno = getGenotypes(obj,varargin)
            p = inputParser;
            p.parse(varargin{:});
            inputs = p.Results;
            
            % get all unique experiment genotypes
            geno = {};
            for i=1:numel(obj.experiment)
                try
                   geno = [geno {obj.experiment(i).geno.label}];
                catch
                end
            end
            geno = unique(geno);
        end
        function [out,data] = heatmap(obj,varargin)
            % inputs
            p = inputParser;
            p.addParameter('type','experiment',@(x) ischar(x) && any(ismember(x,{'experiment','day'})));
            p.addParameter('lights',{'on','off'},@iscellstr);
            p.addParameter('normalized',1,@(x) islogical(x) || isnumeric(x));
            p.addParameter('replaceInf',1,@(x) islogical(x) || isnumeric(x)); % should be take care of in normalization
            p.addParameter('scale',3,@isnumeric);
            p.addParameter('plot',1,@(x) islogical(x) || isnumeric(x));
            p.addParameter('cluster',3,@(x) ismember(x,1:3));
            p.addParameter('save',[],@ischar);
            p.addParameter('include',{},@iscellstr); % combined 'mutant geno'
            p.addParameter('includeDay',{},@iscellstr);
            p.addParameter('includeMutant',{},@iscellstr);
            p.addParameter('includePlate',{},@iscellstr);
            p.addParameter('includeGeno',{},@iscellstr);
            p.addParameter('match',true,@(x) islogical(x) || isnumeric(x));
            p.addParameter('param',{...
                'STIM_mean',...
                'STIM_std',...
                'STIM_propmv',...
                'STIM_nrmMean',... % BASELINEmSTIM_mean
                'POST_mean',...
                'POST_std',...
                'POST_propmv',...
                'STIM1_mean',...
                'STIM1_propmv',...
                'POST1_mean',...
                'POST1_propmv',...
                'STIM2_mean',...
                'STIM2_propmv',...
                'POST2_mean',...
                'POST2_propmv',...
                'STIM3_mean',...
                'STIM3_propmv',...
                'POST3_mean',...
                'POST3_propmv',...
                'STIM4_mean',...
                'STIM4_propmv',...
                'POST4_mean',...
                'POST4_propmv',...
                'STIM5_mean',...
                'STIM5_propmv',...
                'POST5_mean',...
                'POST5_propmv',...
                'STIMmPOST_mean',...
                'S1mS5_mean',...
                },@iscellstr);
            p.parse(varargin{:});
            inputs = p.Results;
            % setup outputs
            out.norm.on = {};
            out.norm.off = {};
            out.rows = {};
            out.cols = {};
            out.data = [];
            if obj.isloaded
                % select data source
                switch inputs.type
                    case 'experiment'
                        data = obj.experiment;
                    case 'day'
                        if ~isempty(obj.day)
                            data = obj.day;
                        else
                            error('Days not combined from experiments');
                        end
                end
                % include plates, mutants, days
                if ~isempty(inputs.includeMutant)
                    data = data(ismember({data.mutant},inputs.includeMutant));
                end
                if ~isempty(inputs.includePlate) && strcmp(inputs.type,'experiment') % doesn't work for data combined across plates
                    data = data(ismember({data.plate},inputs.includePlate));
                end
                if ~isempty(inputs.includeDay)
                    data = data(ismember({data.day},inputs.includeDay)); 
                end
                % organize data
                data = data(ismember({data.lights},inputs.lights));
                for i=1:numel(data)
                    for j=1:numel(data(i).geno)
                        name = sprintf('%s %s',data(i).mutant,data(i).geno(j).label);
                        data(i).name = name;
                        if (isempty(inputs.include) || ismember(name,inputs.include)) && (isempty(inputs.includeGeno) || ismember(data(i).geno(j).label,inputs.includeGeno))
                            l = data(i).lights;
                            % find other indices to make sure rows align
                            switch l
                                case 'on'
                                    lo = 'off';
                                case 'off'
                                    lo = 'on';
                            end
                            idx = find(ismember(out.rows,name));
                            % add to global output
                            if isfield(data(i).geno(j).param,'norm') && ~isempty(data(i).geno(j).param.norm) %&& sum(idx_geno) == 1 && ~isempty(exp(idx_exp).geno(idx_geno).param.norm)
                                % find matched light
                                idx_geno = 0;
                                switch inputs.type
                                    case 'experiment'
                                        idx_data = ismember({data.name},regexprep(data(i).name,l,lo)) & ismember({data.mutant},data(i).mutant);
                                    case 'day'
                                        idx_data = ismember({data.day},data(i).day) & ismember({data.lights},lo) & ismember({data.mutant},data(i).mutant);
                                end
                                if sum(idx_data)
                                    idx_geno = ismember({data(idx_data).geno.label},data(i).geno(j).label);
                                    if sum(idx_geno) == 0
                                       1; 
                                    end
                                end
                                % include current (optionally only include if match found)
                                if ~inputs.match || sum(idx_geno)
                                    if ~isempty(idx)
                                        if numel(out.norm.(l)) < idx
                                            out.norm.(l){idx} = data(i).geno(j).param.norm;
                                        else
                                            out.norm.(l){idx} = [out.norm.(l){idx}; data(i).geno(j).param.norm];
                                        end
                                    else
                                        out.norm.(l){end+1} = data(i).geno(j).param.norm;
                                        out.rows{end+1} = name;
                                        idx = numel(out.rows);
                                    end
                                    if inputs.replaceInf
                                        inf_idx = out.norm.(l){idx} == inf | out.norm.(l){idx} == -inf;
                                        out.norm.(l){idx}(inf_idx) = NaN;
                                    end
                                    if isempty(out.cols)
                                        out.cols = data(i).geno(j).param.columns;
                                    elseif ~all(strcmp(out.cols,data(i).geno(j).param.columns))
                                        error('Column mismatch at %s\t%s',data(i).mutant,data(i).name);
                                    end
                                else
                                    fprintf('SKIPPED: no match found: %s\t%s\n',name);
                                end
                            else
                                fprintf('SKIPPED: norm not calculated or all nan: %s\t%s\n',data(i).mutant,data(i).name);
                            end
                        end
                    end
                end
                % aggregate data
                rm = [];
                for i=1:numel(inputs.lights)
                    l = inputs.lights{i};
                    means.(l) = cellfun(@(x) nanmean(x,1),out.norm.(l),'UniformOutput',0);
                    tmp_rm = cellfun(@(x) all(isnan(x)),means.(l));
                    if isempty(rm)
                        rm = tmp_rm;
                    else
                        rm = rm | tmp_rm;
                    end
                end
                % rm genotypes without data in either
                if sum(rm)
                    fprintf('REMOVING: no data: %s\n',out.rows{rm});
                    out.rows(rm) = [];
                    for i=1:numel(inputs.lights)
                        l = inputs.lights{i};
                        out.norm.(l)(rm,:) = [];
                        means.(l)(rm) = [];
                    end
                end
                if range(structfun(@numel,means))
                    error('Different number of genotypes in ''on'' and ''off''');
                elseif range(struct2array(structfun(@(x) cellfun(@numel,x),means,'UniformOutput',0)))
                    error('Different number of parameters in ''on'' and/or ''off'' data');
                end
                % check param
                if ~isempty(inputs.param)
                    idx = ismember(inputs.param,out.cols);
                    if all(idx)
                        idx = ismember(out.cols,inputs.param);
                        out.cols = out.cols(idx);
                        for i=1:numel(inputs.lights)
                            l = inputs.lights{i};
                            out.norm.(l) = cellfun(@(x) x(:,idx),out.norm.(l),'UniformOutput',0);
                            means.(l) = cellfun(@(x) x(:,idx),means.(l),'UniformOutput',0);
                        end
                    else
                        error('ERROR: %d param missing in loaded data',sum(~idx));
                    end
                end
                % collapse on/off
                out.data = [];
                for i=1:numel(inputs.lights)
                    l = inputs.lights{i};
                    idx = size(out.data,2)+1:size(out.data,2)+numel(out.cols);
                    for j=1:numel(means.(l))
                        out.data(j,idx) = means.(l){j};
                    end
                end
                % column labels light-specific
                cols = out.cols;
                out.cols = {};
                for i=1:numel(inputs.lights)
                    l = inputs.lights{i};
                    out.cols = [out.cols cellfun(@(x) [x '_' l],cols,'UniformOutput',0)];
                end
                % plot
                if inputs.plot
%                     % include
%                     if ~isempty(inputs.includeGeno)
%                         idx = ismember(out.rows,inputs.includeGeno);
%                         out.rows = out.rows(idx);
%                         out.data = out.data(idx,:);
%                         fprintf('%d/%d ''include'' genotypes\n',sum(idx),numel(inputs.includeGeno));
%                     end
                    rm = find(isnan(mean(out.data,1)) | abs(mean(out.data,1)) == inf);
                    for i=1:numel(rm)
                        fprintf('%s removed for inf/nan values\n',out.cols{rm(i)});
                    end
                    out.data(:,rm) = [];
                    out.cols(rm) = [];
                    
                    h = HeatMap(out.data,'columnlabels',regexprep(out.cols,{'_'},' '),'columnlabelsrotate',45,'rowlabels',regexprep(out.rows,'_',' '),'displayrange',inputs.scale);
                    if ~isempty(inputs.save)
                        f = plot(h);
                        saveas(f,[inputs.save '_heatmap']);
                        close(gcf)
                    end
                    h = clustergram(out.data,'columnlabels',regexprep(out.cols,{'_'},' '),'columnlabelsrotate',45,'rowlabels',regexprep(out.rows,'_',' '),'displayrange',inputs.scale,'cluster',inputs.cluster);
                    if ~isempty(inputs.save)
                        f = plot(h);
                        saveas(f,[inputs.save '_clustergram']);
                        close(gcf)
                    end
                end
            else
                error('Experiments not loaded');
            end
        end
    end
    methods (Static)
        function exp = scrapeExp(in,varargin)
            p = inputParser;
            p.addRequired('in',@(x) iscellstr(x) || exist(x,'dir'));
            p.addParameter('match',false,@(x) islogical(x) || isnumeric(x));
            p.addParameter('lights',{'on','off'},@iscellstr);
            p.addParameter('includeMutant',{},@iscellstr);
            p.addParameter('includePlate',{},@iscellstr);
            p.parse(in,varargin{:});
            inputs = p.Results;
            
            % get all .mat 'files' containing animalBehaviorAnalysis objects
            files = {};
            if iscellstr(in)
                if all(cellfun(@(x) exist(x,'file'),in))
                    % make absolute
                    for i=1:numel(in)
                        [fpath,fname,fext] = fileparts(in{i});
                        if isempty(fpath) || ~strcmp(fpath(1),filesep)
                            in{i} = fullfile(pwd,fpath,[fname fext]);
                        end
                    end
                    files = in;
                else
                    error('Not all files inputted exit');
                end
            elseif isdir(in)
                d = dir(in);
                d(~[d.isdir]) = [];
                d(cellfun(@(x) strcmp(x(1),'.'),{d.name})) = [];
                for i=1:numel(d)
                    % get experiments
                    d1 = dir(fullfile(d(i).folder,d(i).name));
                    d1(~[d1.isdir]) = [];
                    d1(cellfun(@(x) strcmp(x(1),'.'),{d1.name})) = [];
                    d1(contains({d1.name},{'combined','reliability'})) = [];
                    for j=1:numel(d1)
                        mat = dir(fullfile(d1(j).folder,d1(j).name,'*.mat'));
                        if numel(mat) == 1
                            files{end+1} = fullfile(mat.folder,mat.name);
                        else
                            fprintf('ERROR: %d .mats\t%s\t%s\n',numel(mat),d(i).name,d1(j).name);
                        end
                    end
                end
                
            end
            % create exp struct from 'files'
            files = unique(files);
            exp = [];
            for i=1:numel(files)
                [fpath,fname] = fileparts(files{i});
                exp(end+1).name = fname;
                [~,exp(end).mutant] = fileparts(fileparts(fpath));
                tmp = strsplit(exp(end).name,'_');
                exp(end).day = tmp{1};
                exp(end).plate = sprintf('%s_%s',tmp{1},tmp{2});
                % include
                if (~isempty(inputs.includeMutant) && ~ismember(exp(end).mutant,inputs.includeMutant)) ||...
                        (~isempty(inputs.includePlate) && ~ismember(exp(end).plate,inputs.includePlate))
                    exp(end) = [];
                    
                else
                    exp(end).mat = files{i};
                    lights = strsplit(regexprep(fname,{'_geno','\.mat'},''),'_');
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
            p.addRequired('geno',@(x) isstruct(x) || iscellstr(x));
            p.addParameter('rename',{},@iscell);
            p.addParameter('suffix',true,@(x) islogical(x) || isnumeric(x));
            p.addParameter('fix',{...
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
            p.addParameter('lower',true,@(x) islogical(x) || isnumeric(x));
            p.parse(geno,varargin{:});
            inputs = p.Results;
            
            % get labels and fix common typos
            if isstruct(geno)
                label = {geno.label};
            else
                label = geno;
            end
            label = cellfun(@strtrim,label,'UniformOutput',0);
            label = regexprep(label,{': ',' :'},':');
            label = regexprep(label,{'/ ',' /'},'/');
            label = regexprep(label,'del ','del');
            
            % check if can be fixed (all genotypes in isolation with one preceded by gene)
            u = unique(label);
            u_new = cell(size(u));
            idx = ismember(u,inputs.fix);
            if sum(idx) == numel(idx) - 1
                % split prefix from genotype
                idx = find(~idx);
                tmp = strsplit(u{idx});
                prefix = lower(join(tmp(1:end-1)));
                prefix = prefix{1};
                for i=1:numel(u)
                    if i==idx
                        u_new{idx} = [prefix ' ' tmp{end}];
                    else
                        u_new{i} = [prefix ' ' u{i}];
                    end
                    % replace in full 'label' list
                    label(strcmp(label,u{i})) = u_new(i);
                end
            end
            if inputs.lower
                label = lower(label);
            end
            label = regexprep(label,'hom','HOM');
            label = regexprep(label,'het','HET');
            label = regexprep(label,'wt','WT');
            % rename manually with string replacement
            if ~isempty(inputs.rename)
                label = startle.renameGeno(label,inputs.rename);
            end
            % only keep suffix
            if inputs.suffix
                for i=1:numel(label)
                    tmp = strsplit(label{i});
                    label{i} = tmp{end};
                end
            end
            if isstruct(geno)
                [geno.label] = label{:};
            else
                geno = label;
            end
        end
        function geno = renameGeno(geno,rename,varargin)
            % inputs
            p = inputParser;
            p.addRequired('geno',@(x) (isstruct(x) && isfield(x,'label')) || iscellstr(x));
            p.addRequired('rename',@iscell);
            p.parse(geno,rename,varargin{:});
            if isstruct(geno)
                label = {geno.label};
            else
                label = geno;
            end
            if ~isempty(rename)
                for j=1:numel(rename)
                    rn = rename{j};
                    if iscell(rn) && numel(rn) == 2
                        label = regexprep(label,rn{1},rn{2});
                    else
                        error('''rename'' should be in {{find1,rep1} {find2,rep2}...{findN,repN}} format');
                    end
                end
            end
            if isstruct(geno)
                [geno.label] = label{:};
            else
                geno = label;
            end
        end
        function geno = normGeno(geno,varargin)
            p = inputParser;
            p.addRequired('geno',@isstruct);
            p.addParameter('normSearch',{'WT' '+'},@iscellstr);
            p.addParameter('normType','z',@(x) ischar(x) && any(ismember(x,{'z','percent'})));
            p.addParameter('propZ',true,@(x) islogical(x) || isnumeric(x));
            p.addParameter('normOutliers',[],@(x) islogical(x) || isnumeric(x));
            p.addParameter('normUpdate',false,@(x) islogical(x) || isnumeric(x));
            p.addParameter('verbose',true,@(x) islogical(x) || isnumeric(x));
            p.parse(geno,varargin{:});
            inputs = p.Results;
            % find normalization index
            name = {geno.label};
            n = zeros(size(name));
            for j=1:numel(name)
                n(j) = max(cellfun(@numel,regexpi(name{j},inputs.normSearch)));
            end
            [max_n,norm_idx] = max(n);
            if ~any(n)
                if inputs.verbose
                    fprintf('\n\t\tNORM FAILED: no matching genotype found');
                end
            elseif sum(max_n == n) > 1
                if inputs.verbose
                    fprintf('\n\t\tNORM FAILED: %d genotypes found',sum(max_n==n));
                end
            else
                % get norm vals
                tmp = geno(norm_idx).param.raw;
                norm.mean = nanmean(tmp,1);
                norm.std = std(tmp,0,1,'omitnan');
                norm.n = size(tmp,1);
                norm.z = (tmp - norm.mean) ./ norm.std;
                % proportion idx
                if inputs.propZ
                    prop_idx = contains(geno(norm_idx).param.columns,'prop');
                else
                    prop_idx = false(size(geno(norm_idx).param.columns)); 
                end
                % exclude norm outliers (don't use prop fields)
                if ~isempty(inputs.normOutliers)
                    rm = abs(norm.z) > inputs.normOutliers;
                    rm(:,prop_idx) = []; 
                    rm = any(rm,2);
                    if any(rm)
                        fprintf('\n\tNORM OUTLIERS: %d/%d excluded',sum(rm),numel(rm));
                        tmp(rm,:) = [];
                        norm.mean = nanmean(tmp,1);
                        norm.std = std(tmp,0,1,'omitnan');
                        if inputs.normUpdate
                            geno(norm_idx).param.raw = tmp;
                            geno(norm_idx).plate(rm) = [];
                            geno(norm_idx).animal_id(rm) = [];
                        end
                    end
                end
                norm.n = size(tmp,1);
                
                % normalize each geno
                invalid = isnan(norm.mean) | norm.mean==inf | norm.mean==-inf;
                if strcmp(inputs.normType,'z')
                    invalid = invalid | isnan(norm.std) | norm.std==0 | norm.std==inf | norm.std==-inf;
                end
                if all(invalid)
                    if inputs.verbose
                        fprintf('\n\t\tNORM FAILED: ');
                        if norm.n == 1 && strcmp(inputs.normType,'z')
                            fprintf('only 1 norm fish');
                        else
                            fprintf('all norm parameters are invalid');
                        end
                    end
                else
                    for j=1:numel(geno)
                        switch inputs.normType
                            case 'z'
                                geno(j).param.norm = (geno(j).param.raw - norm.mean) ./ norm.std;
                                % proportion z
                                if any(prop_idx)
                                    tmp = geno(j).param.raw(:,prop_idx);
                                    tmp(tmp == 1) = 0.99;
                                    tmp(tmp == 0) = 0.01;
                                    tmp_norm = norm.mean(prop_idx);
                                    tmp_norm(tmp_norm == 1) = 0.99;
                                    tmp_norm(tmp_norm == 0) = 0.01;
                                    geno(j).param.norm(:,prop_idx) = (tmp - tmp_norm) ./ sqrt(tmp_norm.*(1-tmp_norm)./norm.n);
                                end
                            case 'percent'
                                geno(j).param.norm = 100.*(geno(j).param.raw - norm.mean) ./ norm.mean;
                        end
                        geno(j).param.data = nanmean(geno(j).param.norm,1);
                        invalid = isnan(geno(j).param.norm) | geno(j).param.norm==inf | geno(j).param.norm==-inf;
                        if all(invalid)
                            if inputs.verbose
                                fprintf('\n\t\tNORM ERROR: %s: all parameters are invalid',geno(j).label);
                            end
                            geno(j).param.norm = [];
                        elseif any(invalid)
                            if inputs.verbose
                                fprintf('\n\t\tNORM WARNING: %s: inf/nan values resulting in %d/%d param',geno(j).label,sum(invalid),numel(invalid));
                            end
                            geno(j).param.norm(invalid) = nan;
                        end
                    end
                end
            end
        end
    end
end
%#ok<*AGROW>
%#ok<*PROPLC>
%#ok<*TRYNC>