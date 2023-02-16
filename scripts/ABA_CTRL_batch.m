function out = ABA_CTRL_batch(data_path,varargin)
% inputs
p = inputParser;
p.addRequired('data_path',@(x) exist(x,'dir'));
p.addParameter('upload',[],@(x) exist(x,'dir'));
p.addParameter('run',true,@(x) islogical(x) || isnumeric(x));
p.addParameter('include',{},@(x) ischar(x) || iscellstr(x));
p.addParameter('overwrite',false,@(x) islogical(x) || isnumeric(x));
p.addParameter('gui',false,@(x) islogical(x) || isnumeric(x));
p.addParameter('stats',true,@(x) islogical(x) || isnumeric(x));
p.addParameter('parallel',true,@(x) islogical(x) || isnumeric(x));
p.addParameter('indiv',true,@(x) islogical(x) || isnumeric(x));
p.addParameter('days',true,@(x) islogical(x) || isnumeric(x));
p.addParameter('combine',true,@(x) islogical(x) || isnumeric(x));
p.addParameter('logFile',[],@(x) ischar(x) || isempty(x));
p.addParameter('liteTs',true,@(x) islogical(x) || isnumeric(x));
p.parse(data_path,varargin{:});
inputs = p.Results;
if ischar(inputs.include)
   inputs.include = {inputs.include};
end
cmd = [];
log = [];

% copy new uploads
if ~isempty(inputs.upload)
   % setup paths
   analysis_path = regexprep(data_path,'data','analysis');
   % go through mutant folders
   d = dir(inputs.upload);
   d(~[d.isdir]) = [];
   d(cellfun(@(x) strcmp(x(1),'.'),{d.name})) = [];
   for i=1:numel(d)
      mut_path = fullfile(d(i).folder,d(i).name);
      d1 = dir(mut_path);
      d1(~[d1.isdir]) = [];
      d1(cellfun(@(x) strcmp(x(1),'.'),{d1.name})) = [];
      % go through experiment folders
      for j=1:numel(d1)
          exp_path = fullfile(d1(j).folder,d1(j).name);
          % delete existing data and analysis
          new_data = fullfile(data_path,d(i).name,d1(j).name);
          new_analysis = fullfile(analysis_path,d(i).name,d1(j).name);
          if exist(new_data,'dir')
             rmdir(new_data,'s');
          end
          if exist(new_analysis,'dir')
             rmdir(new_analysis,'s');
          end
          % copy data
          if ~exist(fileparts(new_data),'dir')
             mkdir(fileparts(new_data));
          end
          copyfile(exp_path,new_data);
      end
   end
   % delete
   try
       rmdir(fullfile(inputs.upload,'*'),'s');
   catch
   end
end

% get top level mutant data dirs
mutant = dir(inputs.data_path);
mutant(~[mutant.isdir]) = [];
mutant(cellfun(@(x) strcmp(x(1),'.'),{mutant.name})) = [];
% scrape indiv exp 'data'
fprintf('Setting up indiv\n');
for i=1:numel(mutant)
    % find each exp folder of mutant_dir
    mutant(i).path = fullfile(data_path,mutant(i).name);
    mutant(i).error = [];
    d = dir(mutant(i).path);
    d(~[d.isdir]) = [];
    d(cellfun(@(x) strcmp(x(1),'.'),{d.name})) = [];
    e = 0;
    exp = [];
    fprintf('\t%s',mutant(i).name);
    for j=1:numel(d)
        % check experiment folder
        r = regexp(d(j).name,'_');
        if numel(r) == 1 || numel(r) == 2
            date = str2double(d(j).name(1:r(1)-1));
            if numel(r) == 2
                plate = d(j).name(r(1)+1:r(2)-1);
            elseif numel(r) == 1
                plate = d(j).name(r(1)+1:end);
            end
            if ~isnan(date) && numel(plate) == 2
                % update 'exp' meta info struct
                fprintf('\n\t\t%s',d(j).name);
                exp_name = d(j).name;
                r = regexp(exp_name,'_');
                if numel(r) == 2
                   exp_name = exp_name(1:r(2)-1);
                end
                exp_path = fullfile(mutant(i).path,d(j).name);
                e = e + 1;
                exp(e).name = exp_name;
                exp(e).path = exp_path;
                exp(e).date = date;
                exp(e).plate = plate;
                exp(e).error = false;
                exp(e).run = true;
                exp(e).msg = {};
                % get/check required files
                files.viewpoint = dir(fullfile(exp_path,[exp_name '*.XLS']));
                if isempty(files.viewpoint)
                   files.viewpoint = dir(fullfile(exp_path,[exp_name '*.xlsx']));
                end
                files.viewpoint(cellfun(@(x) any(~cellfun(@isempty,regexpi(x,{'genotype','timing'}))),{files.viewpoint.name})) = []; % remove other xlsxs
                files.grouping = dir(fullfile(exp_path,[exp_name '*genotype.txt']));
                files.timing = dir(fullfile(exp_path,'*Timing*.txt'));
                files.timing(cellfun(@(x) strcmp(x(1),'.'),{files.timing.name})) = [];
                ff = fields(files);
                for k=1:numel(ff)
                    if numel(files.(ff{k})) == 1
                        exp(e).(ff{k}) = fullfile(files.(ff{k}).folder,files.(ff{k}).name);
                    else
                        % file error
                        exp(e).(ff{k}) = [];
                        exp(e).error = true;
                        msg = sprintf('%d %s files',numel(files.(ff{k})),ff{k});
                        fprintf('\terror: %s',msg);
                        exp(e).msg{end+1} = msg;
                    end
                end
                % load timing info and get lights info (on/off)
                exp(e).lights = [];
                if exist(exp(e).timing,'file')
                    t = animalBehaviorAnalysis.readTimingFile('file',exp(e).timing);
                    if ~t.error
                        try
                            exp(e).lights = t.dataTable.Light{2};
                        catch
                            % timing error
                            exp(e).error = true;
                            msg = 'unexpected timing file structure';
                            fprintf('\terror: %s',msg);
                            exp(e).msg{end+1} = msg;
                        end
                    else
                        % timing error
                       exp(e).error = true;
                       msg = 'couldn''t read timing file';
                       fprintf('\terror: %s',msg);
                       exp(e).msg{end+1} = msg;
                    end
                end
                % add existing analysis .mat
                mat = dir(fullfile(regexprep(mutant(i).path,'data','analysis'),[exp(e).name '_' exp(e).lights],[exp(e).name '_' exp(e).lights '*.mat']));
                mat(contains({mat.name},'_err.mat')) = [];
                mat(contains({mat.name},'_cmd.mat')) = [];
                if numel(mat) == 1
                   exp(e).analysis = fullfile(mat.folder,mat.name);
                   if ~inputs.overwrite
                       fprintf('\tskipped');
                       exp(e).run = false;
                   end
                else
                   exp(e).analysis = [];
                end
            end
        end
    end
    % handle errors
    if ~isempty(exp)
        mutant(i).err = exp([exp.error]==true);
        exp([exp.error]==true) = [];
        mutant(i).exp = exp;
        % move
        for j=1:numel(mutant(i).err)
            [errPath,mutantFolder] = fileparts(mutant(i).path);
            errPath = fullfile(fileparts(errPath),'error',mutantFolder);
            dataPath = mutant(i).err(j).path;
            if ~isdir(errPath)
                mkdir(errPath);
                system(sprintf('chmod 770 "%s"',errPath));
            end
            % copy to error
            code = system(sprintf('rsync -av "%s" "%s" > /dev/null',dataPath,errPath));
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
    end
    % add to indiv cmd batch
    if isempty(inputs.include) || any(strcmp(mutant(i).name,inputs.include))
        % run each experiment seperately
        if inputs.indiv
            for j=1:numel(mutant(i).exp)
                if mutant(i).exp(j).run
                    % setup cmd
                    cmd(end+1).mutant = mutant(i);
                    cmd(end).lights = mutant(i).exp(j).lights;
                    files = [];
                    for k=1:numel(ff)
                        files.(ff{k}) = {mutant(i).exp(j).(ff{k})};
                    end
                    cmd(end).files = files;
                    cmd(end).n = 1;
                    cmd(end).exp = mutant(i).exp(j);
                    cmd(end).savePath = regexprep(mutant(i).path,'data','analysis');
                    cmd(end).projectName = sprintf('%s_%s',mutant(i).exp(j).name,cmd(end).lights);
                end
            end
        end
    else
        fprintf('\n...excluded');
    end
    fprintf('\n\n');
end

% preprocessed
if ~isempty([mutant.err])
    log = [log sprintf('Preprocessing Errors:\n')];
    for i=1:numel(mutant)
        for j=1:numel(mutant(i).err)
            log = [log sprintf('\t%s\t%s',mutant(i).name,mutant(i).err(j).name)];
            for k=1:numel(mutant(i).err(j).msg)
                log = [log sprintf('\t%s',mutant(i).err(j).msg{k})];
            end
            log = [log newline newline];
        end
    end
end

% run all indiv analyses in parallel
if inputs.run && numel(cmd) > 0
    [cmd,tmp_log] = run_cmd(cmd,'parallel',inputs.parallel,'stats',inputs.stats,'liteTs',inputs.liteTs);
    log = [log sprintf('Process individual VSR experiments\t%s\n',datestr(now,'mm/dd/yyyy HH:MM')) tmp_log newline];
    % update mutant struct
    for i=1:numel(cmd)
       midx = ismember({mutant.name},cmd(i).mutant.name);
       eidx = ismember({mutant(midx).exp.path},cmd(i).exp.path);
       % error: remove from mutant.exp and add to mutant.err
       if ~isempty(cmd(i).err) || isempty(cmd(i).exp.analysis)
           mutant(midx).err = [mutant(midx).err cmd(i).exp];
           mutant(midx).exp(eidx) = [];
       else
          mutant(midx).exp(eidx) = cmd(i).exp;
       end
    end
    out.cmd.indiv = cmd;
end

if inputs.combine || inputs.days
    cmd = [];
    % loop through mutant, setup combine cmd
    fprintf('Setting up combination\n')
    for i=1:numel(mutant)
        fprintf('\t%s',mutant(i).name);
        if isempty(inputs.include) || any(ismember(inputs.include,mutant(i).name))
            if inputs.overwrite
                % use all experiments
                days = unique([mutant(i).exp.date]);
            else
                % only use newly run
                days = unique([mutant(i).exp([mutant(i).exp.run] == 1).date]);
            end
            if ~isempty(days)
                % seperately for each light
                lights = unique({mutant(i).exp.lights});
                for j=1:numel(lights)
                    if inputs.days
                        l = lights{j};
                        % each day
                        for k=1:numel(days)
                            exp = mutant(i).exp(ismember({mutant(i).exp.lights},l) & [mutant(i).exp.date] == days(k));
                            if numel(exp) > 1
                                % setup cmd
                                cmd(end+1).mutant = mutant(i);
                                cmd(end).lights = l;
                                cmd(end).files = {exp.analysis};
                                cmd(end).n = numel(exp);
                                cmd(end).exp = exp;
                                cmd(end).savePath = regexprep(mutant(i).path,'data','analysis');
                                cmd(end).projectName = sprintf('%d_%s',days(k),l);
                                fprintf('\n\t\t%s',cmd(end).projectName);
                            end
                        end
                    end
                    % all combined
                    if inputs.combine && (inputs.overwrite || sum([mutant(i).exp.run])) && any(ismember({mutant(i).exp.lights},l))
                        % add all cmd
                        exp = mutant(i).exp(ismember({mutant(i).exp.lights},l));
                        if numel(exp) > 1
                            cmd(end+1).mutant = mutant(i);
                            cmd(end).lights = l;
                            cmd(end).files = {exp.analysis};
                            cmd(end).n = numel(numel(exp));
                            cmd(end).exp = exp;
                            cmd(end).savePath = regexprep(mutant(i).path,'data','analysis');
                            cmd(end).projectName = sprintf('combined_%s',l);
                            fprintf('\n\t\t%s',cmd(end).projectName);
                        end
                    end
                end
            else
               fprintf('\tno new');
            end
        else
            fprintf('\texcluded');
       end
       fprintf('\n');
    end
   if ~isempty(cmd)
       [cmd,tmp_log] = run_cmd(cmd,'parallel',inputs.parallel,'stats',inputs.stats,'liteTs',inputs.liteTs);
       log = [log sprintf('Combine VSR experiments\t%s\n',datestr(now,'mm/dd/yyyy HH:MM')) tmp_log];
       out.cmd.combined = cmd;
   end
end

% log / output
out.mutant = mutant;
out.log = log;
if ~isempty(inputs.logFile)
    [logPath,logName] = fileparts(inputs.logFile);
    if ~isempty(out.log)
        fid = fopen(fullfile(logPath,[logName '.log']),'w');
        fprintf(fid,'%s\n',log);
        fclose(fid);
    end
    save(fullfile(logPath,logName),'out');
end
end

function [cmd,log] = run_cmd(cmd,varargin)
p = inputParser;
p.addRequired('cmd',@isstruct);
p.addParameter('parallel',true,@(x) islogical(x) || isnumeric(x));
p.addParameter('stats',true,@(x) islogical(x) || isnumeric(x));
p.addParameter('liteTs',true,@(x) islogical(x) || isnumeric(x));
p.addParameter('gui',false,@(x) islogical(x) || isnumeric(x));
p.parse(cmd,varargin{:});
inputs = p.Results;
% run_aba in parallel
if inputs.parallel
    start_parpool;
    clear f
    for i=1:numel(cmd)
        cmd(i).err = [];
        f(i) = parfeval(@run_aba,1,cmd(i),'stats',inputs.stats,'liteTs',inputs.liteTs);
    end
    % fetch and waitbar
    if inputs.gui
        h = waitbar(0,sprintf('Analyzing %d datasets',numel(cmd)));
        count = 0;
    end
    for i=1:numel(cmd)
        try
            [idx,val] = fetchNext(f);
            cmd(idx) = val;
        catch
        end
        if inputs.gui
            count = count + 1;
            waitbar(count/numel(subj),h);
        end
    end
    % close waitbar
    if inputs.gui
        try
            close(h);
        catch
            delete(h);
        end
    end
    % get diary
    for i=1:numel(cmd)
        str = f(i).Diary;
        str = strsplit(str,'\n');
        str(cellfun(@isempty,str)) = [];
        cmd(i).diary = str;
%         for j=1:numel(str)
%             fprintf('%s\n',str{j});
%         end
%         fprintf('\n');
    end
else
    for i=1:numel(cmd)
        cmd(i).err = [];
        cmd(i) = run_aba(cmd(i),'stats',inputs.stats,'liteTs',inputs.liteTs);
    end
end

% processed
log = sprintf('Processed:\n');
for i=1:numel(cmd)
    log = [log sprintf('\t%s\t%s',cmd(i).mutant.name,cmd(i).projectName)];
    if ~isempty(cmd(i).err)
        log = [log sprintf('\tERROR')];
    end
    log = [log newline];
end
end

%#ok<*AGROW>
