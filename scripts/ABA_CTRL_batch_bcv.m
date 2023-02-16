function [cmd,mutant] = ABA_CTRL_batch_bcv(data_path,varargin)
% inputs
p = inputParser;
p.addRequired('data_path',@(x) exist(x,'dir'));
p.addParameter('upload',[],@(x) exist(x,'dir'));
p.addParameter('run',true,@(x) islogical(x) || isnumeric(x));
p.addParameter('include',{},@iscellstr);
p.addParameter('overwrite',false,@(x) islogical(x) || isnumeric(x));
p.addParameter('gui',false,@(x) islogical(x) || isnumeric(x));
p.addParameter('stats',true,@(x) islogical(x) || isnumeric(x));
p.addParameter('parallel',true,@(x) islogical(x) || isnumeric(x));
p.addParameter('combine',true,@(x) islogical(x) || isnumeric(x));
p.addParameter('indiv',true,@(x) islogical(x) || isnumeric(x));
p.addParameter('summary',[],@(x) ischar(x) || isempty(x));
p.addParameter('liteTs',true,@(x) islogical(x) || isnumeric(x));
p.parse(data_path,varargin{:});
inputs = p.Results;

cmd = [];

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
   unix(sprintf('rm -r %s/* 2> /dev/null',inputs.upload));
end

% get top level mutant dirs
mutant = dir(inputs.data_path);
mutant(~[mutant.isdir]) = [];
mutant(cellfun(@(x) strcmp(x(1),'.'),{mutant.name})) = [];
for i=1:numel(mutant)
    % find exp folders of mutant_dir
    mutant(i).path = fullfile(data_path,mutant(i).name);
    mutant(i).error = [];
    d = dir(mutant(i).path);
    d(~[d.isdir]) = [];
    d(cellfun(@(x) strcmp(x(1),'.'),{d.name})) = [];
    e = 0;
    exp = [];
    fprintf('%s',mutant(i).name);
    for j=1:numel(d)
        % check experiment folder naming
        r = regexp(d(j).name,'_');
        if numel(r) == 1 || numel(r) == 2
            date = str2double(d(j).name(1:r(1)-1));
            if numel(r) == 2
                plate = d(j).name(r(1)+1:r(2)-1);
            elseif numel(r) == 1
                plate = d(j).name(r(1)+1:end);
            end
            if ~isnan(date) && numel(plate) == 2
                % update 'ex'p meta info struct
                fprintf('\n\t%s',d(j).name);
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
                    try
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
                    catch
                       1; 
                    end
                end
                % check existing
                mat = dir(fullfile(regexprep(mutant(i).path,'data','analysis'),[exp(e).name '_' exp(e).lights],[exp(e).name '_' exp(e).lights '*.mat']));
                mat(contains({mat.name},'_err.mat')) = [];
                if ~isempty(mat) && ~inputs.overwrite
                    fprintf('\tskipped');
                    exp(e).run = false;
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
            % copy
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
    
    % add to cmd batch
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
        
        % combine experiments
        if inputs.combine && ~isempty(mutant(i).exp)
            % combine lights "on" and "off" seperately
            lights = {mutant(i).exp.lights};
            lights(cellfun(@isempty,lights)) = [];
            lights = unique(lights);
            for j=1:numel(lights)
                idx = strcmp({mutant(i).exp.lights},lights{j});
                idxRun = idx & [mutant(i).exp.run];
                projectName = sprintf('combined_%s',lights{j});
                % check existing
                mat = fullfile(regexprep(mutant(i).path,'data','analysis'),projectName,[projectName '_err.mat']);
                if (exist(mat,'file') || inputs.overwrite || sum(idxRun)) && sum(idx) > 1 
                    for k=1:numel(ff)
                        files.(ff{k}) = {mutant(i).exp(idx).(ff{k})};
                    end
                    cmd(end+1).mutant = mutant(i);
                    cmd(end).lights = lights{j};
                    cmd(end).files = files;
                    cmd(end).n = sum(idx);
                    cmd(end).exp = mutant(i).exp(idx);
                    cmd(end).savePath = regexprep(mutant(i).path,'data','analysis');
                    cmd(end).projectName = projectName;
                    fprintf('\n\t%s',cmd(end).projectName);
                else
                    fprintf('\n\t%s\tskipped',projectName);
                end
            end
        end
    else
        fprintf('\n...excluded');
    end
    fprintf('\n\n');
end

% run all animalBehavior analyses in parallel
if inputs.run && numel(cmd) > 0
    % run fcn in parallel
    if inputs.parallel
        start_parpool;
        clear f
        for i=1:numel(cmd)
            cmd(i).err = [];
            f(i) = parfeval(@run_aba,1,cmd(i),'stats',inputs.stats,'liteTs',inputs.liteTs);
        end
        % fetch and waitbar
        if inputs.gui
            h = waitbar(0,sprintf('Analyzing %d combined data sets',numel(cmd)));
            count = 0;
        end
        for i=1:numel(cmd)
            try
                [idx,val] = fetchNext(f);
                cmd(idx) = val;
            catch
                1;
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
        % print diary
        for i=1:numel(cmd)
            str = f(i).Diary;
            str = strsplit(str,'\n');
            str(cellfun(@isempty,str)) = [];
            for j=1:numel(str)
                fprintf('%s\n',str{j});
            end
            fprintf('\n');
        end
    else
        for i=1:numel(cmd)
            cmd(i).err = [];
            cmd(i) = run_aba(cmd(i),'stats',inputs.stats,'liteTs',inputs.liteTs);
        end
    end
    % summary output
    if ~isempty(inputs.summary)
       fid = fopen(inputs.summary,'w');
       fprintf(fid,'New VSR experiments processed\t%s\n\n',datestr(now,'mm/dd/yyyy HH:MM'));
       % preprocessed
       if ~isempty([mutant.err])
           fprintf(fid,'Preprocessing Errors:\n');
           for i=1:numel(mutant)
              for j=1:numel(mutant(i).err)
                 fprintf(fid,'\t%s\t%s',mutant(i).name,mutant(i).err(j).name);
                 for k=1:numel(mutant(i).err(j).msg)
                    fprintf(fid,'\t%s',mutant(i).err(j).msg{k}); 
                 end
                 fprintf(fid,'\n');
              end
           end
       end
       % processed
       fprintf(fid,'Processed:\n');
       for i=1:numel(cmd)
           fprintf(fid,'\t%s\t%s',cmd(i).mutant.name,cmd(i).projectName);
           if ~isempty(cmd(i).err)
              fprintf(fid,'\tERROR'); 
           end
           fprintf(fid,'\n');
       end
       fclose(fid);
    end
end
end

%#ok<*AGROW>