% script to combine multiple animalBehavioralAnalysis objects, output
% timeseries plots, reliability (z-score relative to WT) csv, and heatmap
%
% SurveyBott, 2019, info@surveybott.com
%#ok<*AGROW>
function data = combinedReliability(in,varargin)
wt = fullfile(fileparts(mfilename('fullpath')),'wt.mat');
% inputs
p = inputParser;
p.addRequired('in',@isdir);
p.addParameter('out',fullfile(in,'reliability'),@ischar);
p.addParameter('ts',true,@(x) islogical(x) || isnumeric(x));
p.addParameter('zscore',true,@(x) islogical(x) || isnumeric(x));
p.addParameter('mat',true,@(x) islogical(x) || isnumeric(x));
p.addParameter('heatmap',true,@(x) islogical(x) || isnumeric(x));
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
p.addParameter('wt',wt,@(x) exist(x,'file'));
p.addParameter('addNorm',true,@(x) islogical(x) || isnumeric(x));
p.parse(in,varargin{:});
inputs = p.Results;

% load WT norm data, compute z-score values
wt = load(inputs.wt);
lights = {'on','off'};
for i=1:numel(lights)
    if all(ismember(inputs.param,wt.uni.(lights{i}).Properties.VariableNames))
       columns = ismember(wt.uni.(lights{i}).Properties.VariableNames,inputs.param);
       mat = wt.uni.(lights{i}){:,columns};
       wt.param.(lights{i}).mean = nanmean(mat,1);
       wt.param.(lights{i}).std = std(mat,0,1,'omitnan');
       wt.param.(lights{i}).cols = wt.uni.(lights{i}).Properties.VariableNames(columns);
    else
       error('Not all parameters found in wt.uni.%s',lights{i}); 
    end
end

geno = {'wt','hom'};

% run for both lights on and off
for i=1:numel(lights)
    for g=1:numel(geno)
        data.(lights{i}).(geno{g}).mean = [];
        data.(lights{i}).(geno{g}).std = [];
        data.(lights{i}).(geno{g}).z = [];
        data.(lights{i}).(geno{g}).rows = {};
        data.(lights{i}).(geno{g}).cols = wt.param.(lights{i}).cols;
    end
   
   l = lights{i};
   % find .mat
   mat = dir(fullfile(in,sprintf('*_*_%s',l),sprintf('*_*_%s.mat',l)));
   % add experiment days
   for j=1:numel(mat)
       day = strsplit(mat(j).name,'_');
       mat(j).day = day{1};
   end
   % combine all experiments on the same day
   day = unique({mat.day});
   for j=1:numel(day)
      % load aba objects and combine obj.statistics.ts and obj.statistics.uni
      dayMat = mat(strcmp({mat.day},day{j}));
      for g=1:numel(geno)
         uni.(geno{g}) = [];
         ts.(geno{g}) = [];
      end
      for k=1:numel(dayMat)
          tmp = load(fullfile(dayMat(k).folder,dayMat(k).name));
          if isfield(tmp,'n') && isa(tmp.n,'animalBehaviorAnalysis')
              obj = tmp.n;
              % pull HOM & WT
              id = idGeno(unique(obj.data.group_data.group_label2));
              for g=1:numel(geno)
                  if ~isempty(id.(geno{g}))
                      uni.(geno{g}) = [uni.(geno{g}); obj.statistics.uni(strcmp(obj.statistics.uni.group_label2,id.(geno{g})),:)];
                      ts.(geno{g}) = [ts.(geno{g}); obj.statistics.ts(strcmp(obj.statistics.ts.group_label2,id.(geno{g})),:)];
                  else
                      fprintf('%s: no %s found!\n',dayMat(k).name,upper(geno{g}));
                  end
              end
          end
      end
      % output WT timeseries
      if ~isempty(ts.wt)
          ts.wt.group_label2(:) = {'WT'};
          if inputs.ts && ~isempty(inputs.out)
              obj.meta.save_location = fullfile(inputs.out,sprintf('%s_%s',day{j},l));
              if ~isdir(obj.meta.save_location)
                  mkdir(obj.meta.save_location);
              end
              if inputs.addNorm
                  ts.wt = [ts.wt; wt.ts.(lights{i})];
              end
              obj.statistics.ts = ts.wt;
              try
                  obj.basicTimeSeriesPlots;
              catch err
                  fprintf('%s_%s: timeseries plot error\n',day{j},lights{i});
              end
          end
      end
      % compute parameter mean and z-scores (relative to WT norm)
      for g=1:numel(geno)
          if ~isempty(uni.(geno{g}))
              % pull columns in the right order to match WT
              cols = data.(lights{i}).(geno{g}).cols;
              idx = zeros(size(cols));
              for k=1:numel(cols)
                  tmp = find(strcmp(uni.(geno{g}).Properties.VariableNames,cols{k}));
                  if numel(tmp) == 1
                      idx(k) = tmp;
                  else
                      error('%s_%s (%s): %s column not found in obj.statistics.uni table',day{j},lights{i},geno{g},cols{k});
                  end
              end
              % calculate values, normalize
              tmp = uni.(geno{g}){:,idx};
              data.(lights{i}).(geno{g}).mean(end+1,:) = nanmean(tmp,1);
              data.(lights{i}).(geno{g}).std(end+1,:) = std(tmp,0,1,'omitnan');
              data.(lights{i}).(geno{g}).rows{end+1} = day{j};
              switch geno{g}
                  case 'wt' % normalize by WT norm
                      data.(lights{i}).(geno{g}).z(end+1,:) = nanmean((tmp - wt.param.(lights{i}).mean) ./ wt.param.(lights{i}).std,1);
                  case 'hom' % normalize by WT combined on experiment day
                      wt_idx = strcmp(data.(lights{i}).wt.rows,day{j});
                      if sum(wt_idx) == 1
                          data.(lights{i}).(geno{g}).z(end+1,:) = nanmean((tmp - data.(lights{i}).wt.mean(wt_idx,:)) ./ data.(lights{i}).wt.std(wt_idx,:),1);
                      else
                          data.(lights{i}).(geno{g}).z(end+1,:) = nan(1,size(tmp,2));
                      end
              end
          end
      end
   end
   % output heatmap
   for g=1:numel(geno)
       if ~isempty(uni.(geno{g}))
           if ~isempty(data.(lights{i}).(geno{g}).z) && ~isempty(inputs.out) && inputs.heatmap
               warning('off','MATLAB:Figure:FigureSavedToMATFile');
               switch geno{g}
                   case 'wt'
                       title = sprintf('WT Z-Scores Relative to WT Norm (%s)',lights{i});
                   case 'hom'
                       title = sprintf('HOM Z-Scores Relative to WT (%s)',lights{i});
               end
               hm = HeatMap(data.(lights{i}).(geno{g}).z,'ColumnLabels',regexprep(data.(lights{i}).(geno{g}).cols,'_',' '),'ColumnLabelsRotate',45,'RowLabels',regexprep(data.(lights{i}).(geno{g}).rows,'_',' '),'DisplayRange',3,'title',title);
               f = figure;
               ax = hm.plot(f);
               colorbar('peer',ax);
               if ~isdir(inputs.out)
                   mkdir(inputs.out);
               end
               saveas(ax,fullfile(inputs.out,sprintf('heatmap_%s_%s.fig',geno{g},lights{i})));
               close all hidden
           end
           % output z-scores
           if ~isempty(data.(lights{i}).(geno{g}).z) && ~isempty(inputs.out) && inputs.zscore
               if ~isdir(inputs.out)
                   mkdir(inputs.out);
               end
               fid = fopen(fullfile(inputs.out,sprintf('z_avg_%s_%s.tsv',geno{g},lights{i})),'w');
               for j=1:size(data.(lights{i}).(geno{g}).z,1)
                   fprintf(fid,'%s_%s\t%.3f\n',data.(lights{i}).(geno{g}).rows{j},lights{i},nanmean(abs(data.(lights{i}).(geno{g}).z(j,:))));
               end
               fclose(fid);
           end
       end
   end
end

if ~isempty(inputs.out) && inputs.mat
    if isdir(inputs.out) % don't output if no other outputs
        save(fullfile(inputs.out,'reliability'),'data');
    end
end
end