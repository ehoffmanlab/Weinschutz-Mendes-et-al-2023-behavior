% out = ABAERRREPORT collates and reports errors errors from the
% animalBehaviorAnalysis pipeline (ABA_CTRL_batch)
% 
% Required
%   'path_in' location of _err.mat aba files in zebrafish startle/analysis directory or same dir structure.
%       default: pwd       
%
% Parameters
%   'depth' depth of *_err.mat files relative to 'path_in'
%       default: 2
%
% Output
%  'out.mat'    - _err.mat files and locations (from dir function)
%  'out.report' - outputs from MException.getReports
%  'out.unique' - unique errors in 'out.report'
%  'out.idx'    - indices in 'out.mat' & 'out.report' for each 'out.unique'
%
% SurveyBott, 2020, info@surveybott.com

function out = abaErrReport(path_in,varargin)
p = inputParser;
p.addRequired('path_in',@isfolder);
p.addParameter('depth',false);
if ~exist('path_in','var')
   path_in = pwd; 
end
p.parse(path_in,varargin{:});
inputs = p.Results;

% get error mats
str = '*_err.mat';
for i=1:inputs.depth
    str = ['*/' str];
end
out.mat=dir(str);
out.report = cell(size(out.mat));
out.unique = {};
out.idx = cell(size(out.mat));
fprintf('  ');
for i=1:numel(out.mat)
    msg = sprintf('%d/%d...',i,numel(out.mat));
    fprintf('%s',msg);
    del = repmat(sprintf('\b'),1,numel(msg)-3);
    try
        cmd = [];
        load(fullfile(out.mat(i).folder,out.mat(i).name));
        report{i} = cmd.err.getReport;
        fprintf('\b\b\b');
        pause(0.1);
        fprintf('%s',del);
    catch
        fprintf(' error getting error lol\n');
        pause(1);
    end
end

% get and print unique
err = unique(report);
dash = sprintf('----------------------------------------\n');
fprintf('\b%s%d Unique Errors in %d\n%s\n',dash,numel(err),numel(report),dash);
for i=1:numel(err)
    out.idx{i} = find(ismember(report,err{i}));
    fprintf('%s%d\t''report'' indices: ',dash,i);
    fprintf('%d ',out.idx{i});
    fprintf('\n');
    fprintf('%s%s\n\n',dash,err{i})
end
end