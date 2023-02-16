function out = fixTrailing(in,varargin)
p = inputParser;
p.addRequired('in',@(x) exist(x,'file'));
p.addParameter('delimiter','\t');
p.addParameter('out','');
p.parse(in,varargin{:});
inputs = p.Results;
% setup output
if isempty(inputs.out)
    [fpath,fname,fext] = fileparts(in);
    inputs.out = fullfile(fpath,[fname '_fix' fext]);
end

% read input
txt = {};
fid = fopen(in,'r');
while ~feof(fid)
    txt{end+1} = fgetl(fid);
end
fclose(fid);

% fix trailing and output
fid = fopen(inputs.out,'w');
n = 0;
for i=1:numel(txt)
    s = strsplit(txt{i},inputs.delimiter,'CollapseDelimiters',false);
    rm = false(size(s));
    % remove trailing after end of header
    for j=numel(s):-1:1
        if isempty(s{j}) && j > n
           rm(j) = true; 
        else
           break 
        end
    end
    s(rm) = [];
    if i==1
       n = numel(s);
    end
    for j=1:numel(s)
       fprintf(fid,'%s',s{j});
       if isempty(s{j})
          %fprintf(fid,inputs.delimiter); 
       end
       if j<numel(s)
          fprintf(fid,inputs.delimiter);
       else
          fprintf(fid,'\n'); 
       end
    end
end
fclose(fid);
out = inputs.out;
end