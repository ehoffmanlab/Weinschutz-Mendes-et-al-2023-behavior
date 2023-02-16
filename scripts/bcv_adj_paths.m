currdir = fullfile(pwd)

cmd.savePath = currdir 

[PATHSTR,NAME,EXT] = fileparts(cmd.files.viewpoint{1})
cmd.files.viewpoint{1} = fullfile(pwd,[NAME EXT])
[PATHSTR,NAME,EXT] = fileparts(cmd.files.grouping{1})
cmd.files.grouping{1} = fullfile(pwd,[NAME EXT])
[PATHSTR,NAME,EXT] = fileparts(cmd.files.timing{1})
cmd.files.timing{1} = fullfile(pwd,[NAME EXT])

%% merge version

currdir = fullfile(pwd, '\..')

cmd.savePath = currdir 
for i = 1 :  numel(cmd.files)
    [PATHSTR,NAME,EXT] = fileparts(cmd.files{i})
    cmd.files{i} = fullfile(pwd,[NAME EXT])
end

%%
cmd.savePath =  'C:\Users\bv8\Documents\Projects\hoffman\debug_0423\' ;

cmd.files{1} =    'C:\Users\bv8\Documents\Projects\hoffman\debug_0423\/190626_0A_on/190626_0A_on.mat'
cmd.files{2} =   'C:\Users\bv8\Documents\Projects\hoffman\debug_0423\190626_0B_on/190626_0B_on.mat'
cmd.files{3} =   'C:\Users\bv8\Documents\Projects\hoffman\debug_0423\190626_0C_on/190626_0C_on.mat'
cmd.files{4} =   'C:\Users\bv8\Documents\Projects\hoffman\debug_0423\190626_0D_on/190626_0D_on.mat'


