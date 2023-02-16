function xls2txt(in,varargin)
p = inputParser;
p.addRequired('in',@(x) exist(x,'file'));
p.addParameter('delimiter','\t');
p.addParameter('out','',@ischar);
p.parse(in,varargin{:});
inputs = p.Results;

% load input
[~,~,txt] = xlsread(in);

% setup output
[fpath,fname] = fileparts(in);
if isempty(inputs.out)
   inputs.out = fullfile(fpath,[fname '.txt']); 
end
fid = fopen(inputs.out,'w');

% write
for i=1:size(txt,1)
   for j=1:size(txt,2)
      val = txt{i,j};
      if isnumeric(val)
         if round(val) == val
            val = sprintf('%d',val);
         else
            val = sprintf('%f',val); 
         end
      end
      fprintf(fid,val);
      if j < size(txt,2)
         fprintf(fid,inputs.delimiter);
      else
         fprintf(fid,'\n'); 
      end
   end
end
fclose(fid);

end