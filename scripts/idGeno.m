% function to identify WT and HOM genotypes from a list
function out = idGeno(in,varargin)
p = inputParser;
p.addRequired('in',@iscellstr);
p.parse(in,varargin{:});
inputs = p.Results;

out.wt = [];
out.hom = [];

% wt
idx = cellfun(@numel,regexpi(in,'wt')) + cellfun(@numel,regexpi(in,'+'));
idx = find(max(idx) == idx);
if numel(idx) == 1
    out.wt = in{idx};
end

% hom
% first try to id most 'HOM'
idx = cellfun(@numel,regexpi(in,'hom'));
idx = find(max(idx) == idx);
if numel(idx) == 1
    out.hom = in{idx};
else
    % then search for only geno without '+'
    idx = cellfun(@numel,regexpi(in,'+'));
    if sum(idx == 0) == 1 && sum(idx > 0) == numel(idx)-1
        out.hom = in{idx==0};
    end
end
end