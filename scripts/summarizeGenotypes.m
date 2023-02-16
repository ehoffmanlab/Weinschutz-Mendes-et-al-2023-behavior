% Pull group info from multiple grouping files. Uses unix 'find' to get all *genotype.txt files 
%
% Jeff Eilbott, 2017, jeilbott@surveybott.com
function [genotypes,G] = summarizeGenotypes(path,varargin)
p = inputParser;
p.addRequired('path',@isdir);
p.addParameter('findStr','*genotype.txt',@ischar);
p.addParameter('tableCol','group_label2',@ischar);
p.addParameter('print',1,@(x) islogical(x) || isboolean(x));
p.parse(path,varargin{:});
inputs = p.Results;
if isempty(inputs.path)
   error('''path'' parameter is required'); 
end

% find files
[~,files] = unix(sprintf('find "%s" -name "%s"',inputs.path,inputs.findStr));
files = strsplit(files,'\n');
files(cellfun(@isempty,files)) = [];

% parse files, get unique genotypes and count individuals in them
geno_u = {};
genotypes = [];
for i=1:numel(files)
    G(i,1).file = files{i};
    g = animalBehaviorAnalysis.readGroupingFile(files{i});
    G(i).dataTable = g.dataTable;
    geno_u = unique(G(i).dataTable.(inputs.tableCol));
    for j=1:numel(geno_u)
       count = sum(strcmp(G(i).dataTable.(inputs.tableCol),geno_u{j}));
       if ~isstruct(genotypes) || ~isfield(genotypes,'name') || ~any(strcmp({genotypes.name},geno_u{j}))
          genotypes(numel(genotypes)+1,1).name = geno_u{j};
          genotypes(end).count = count;
       else
           idx = strcmp({genotypes.name},geno_u{j});
           genotypes(idx).count = genotypes(idx).count + count;
       end
    end
end

% print info
if inputs.print
   for i=1:numel(genotypes)
      fprintf('%d\t%s\t%d\n',i,genotypes(i).name,genotypes(i).count);
   end
end
end