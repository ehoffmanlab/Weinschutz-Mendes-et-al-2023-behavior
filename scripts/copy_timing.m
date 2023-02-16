% COPY_TIMING(mutants,cmd,...) copies template timing files into data directories
%
% Inputs
%   REQUIRED
%   'mutants'   struct of mutant directories with .exp fields
%   'cmd'      struct that informs what timing files to copy where
%       fields:
%           'mutant'    maps to mutants(i).name
%           'exp'       maps to mutants(i).exp(j).name
%           'lights'    {'on','off'}
%   OPTIONAL PARAMETERS
%   'timing_dir'    location where txt template files are
%   'overwire'      0|1|true|false - overwrite existing timing txt?
%   'verbose'       0|1|true|false
%       -
%
% Jeff Eilbott, 2017, jeilbott@surveybott.com

function mutants = copy_timing(mutants,cmd,varargin)
% inputs
lights = {'on','off'};
p = inputParser;
p.addRequired('mutants',@(x) isstruct(x) && all(isfield(x,{'name','exp'})) && ...
    all(isfield(x(1).exp,{'name','path'})));
p.addRequired('cmd',@(x) isstruct(x) && all(isfield(x,{'mutant','exp','lights'})) && ...
    all(cellfun(@(y) any(strcmp(y,lights)),{x.lights})));
p.addParameter('timing_dir','/ysm-gpfs/project/ejh22/bioinfo/startle/info/timing_files',@isdir);
p.addParameter('overwrite',false,@(x) islogical(x) || isnumeric(x));;
p.addParameter('verbose',true,@(x) islogical(x) || isnumeric(x));
p.parse(mutants,cmd,varargin{:});
inputs = p.Results;

% set up timing files for lights 'on' and 'off'
for i=1:numel(lights)
   timing.(lights{i}) = fullfile(inputs.timing_dir,sprintf('lights%s_Timing.txt',lights{i})); 
end

% loop over copy cmds 
for i=1:numel(cmd)
    if inputs.verbose
        fprintf('(%d/%d)\t%s\t%s\t%s\t',i,numel(cmd),cmd(i).mutant,cmd(i).exp,cmd(i).lights);
    end
    m_idx = strcmp(cmd(i).mutant,{mutants.name});
    if sum(m_idx) == 1
        e_idx = strcmp(cmd(i).exp,{mutants(m_idx).exp.name});
        if sum(e_idx) == 1 
            timing_file = fullfile(mutants(m_idx).exp(e_idx).path,sprintf('%sTiming.txt',cmd(i).exp));
            if exist(timing_file,'file') && ~inputs.overwrite && inputs.verbose
                fprintf('EXISTS and overwrite isn''t set');
            else
                copyfile(timing.(cmd(i).lights),timing_file);
                mutants(m_idx).exp(e_idx).timing = timing_file;
            end
        elseif inputs.verbose
            fprintf('ERROR: %d exp dirs found',sum(e_idx));
        end
    elseif inputs.verbose
        fprintf('ERROR: %d mutant dirs found',sum(m_idx));
    end
    if inputs.verbose
       fprintf('\n'); 
    end
end

end