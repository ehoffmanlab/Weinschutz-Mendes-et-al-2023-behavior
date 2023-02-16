function liteTimeSeries(filename,varargin)
% inputs
p = inputParser;
p.addRequired('filename',@(x) exist(x,'file'));
p.addParameter('keep',{'WT:WT','HET:HET','HOM:HOM'},@iscellstr);
p.addParameter('save','',@ischar);
p.parse(filename,varargin{:});
inputs = p.Results;

% open fig, get legend, and names
f = open(filename);
ax = gca;
names = ax.Legend.String;
idx = contains(names,inputs.keep);
rm_idx = find(~idx);

if sum(idx) == numel(inputs.keep)
    % remove other lines and faces (two axis children / group)
    rm = false(size(ax.Children));
    for i=1:numel(rm_idx)
        rm([rm_idx(i)*2-1 rm_idx(i)*2]) = true;
    end
    delete(ax.Children(rm));
    
    % rearrange for consistency
    names = names(idx);
    idx = zeros(1,numel(inputs.keep)*2);
    for i=1:numel(inputs.keep)
       tmp = find(contains(names,inputs.keep{i}));
       idx(i*2-1) = tmp*2-1;
       idx(i*2) = tmp*2;
    end
    ax.Children = ax.Children(idx);
    % change colors for consistency
    c = lines;
    for i=1:numel(inputs.keep)
       ax.Children(i*2-1).Color = c(i,:);
       ax.Children(i*2).FaceColor = c(i,:);
    end
    
    % update legend
    legend(ax.Children(1:2:end),inputs.keep);
    
    % save
    if ~isempty(inputs.save)
        saveFile = inputs.save;
    else
        [fpath,fname,fext] = fileparts(filename);
        saveFile = fullfile(fpath,[fname '_lite' fext]);
    end
    [fpath,fname,fext] = fileparts(saveFile);
    if isempty(fext)
        saveFile = fullfile(fpath,[fname '.fig']);
    end
    saveas(f,saveFile);
end
close(f);

end