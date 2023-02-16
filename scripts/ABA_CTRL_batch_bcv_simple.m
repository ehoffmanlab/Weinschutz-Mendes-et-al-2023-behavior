cd 'Z:\721430-MIMED_Geriatrics-DMIC06\Biostat\vanderwyk\projects\hoffman\data' 

stats = 1 
liteTs = 0

cmd = load('cmd.mat')
cmd = cmd.cmd
ABA_CTRL_batch_bcv('Z:\721430-MIMED_Geriatrics-DMIC06\Biostat\vanderwyk\projects\hoffman\data',...
                'run',true,...
                'include',{},...
                'overwrite',true,...
                'gui',false,...
                'stats',true,...
                'parallel',false,...
                'combine',false,...
                'indiv',true,...
                'summary',[],...
                'liteTs',false,...
                'reliability',false) ;
% raba = run_aba(cmd,'stats',stats,'liteTs',liteTs)


%% testing merge
projectName = 'merge';
savepath = 'Z:\721430-MIMED_Geriatrics-DMIC06\Biostat\vanderwyk\projects\hoffman\analysis\chd8_del5\';
mergelist = {'Z:\721430-MIMED_Geriatrics-DMIC06\Biostat\vanderwyk\projects\hoffman\analysis\chd8_del5\180702_0A_off\180702_0A_off.mat'
    'Z:\721430-MIMED_Geriatrics-DMIC06\Biostat\vanderwyk\projects\hoffman\analysis\chd8_del5\180702_0B_off\180702_0B_off.mat'
    'Z:\721430-MIMED_Geriatrics-DMIC06\Biostat\vanderwyk\projects\hoffman\analysis\chd8_del5\180716_0C_off\180716_0C_off.mat'
    'Z:\721430-MIMED_Geriatrics-DMIC06\Biostat\vanderwyk\projects\hoffman\analysis\chd8_del5\180813_0A_off\180813_0A_off.mat'
    'Z:\721430-MIMED_Geriatrics-DMIC06\Biostat\vanderwyk\projects\hoffman\analysis\chd8_del5\180813_0B_off\180813_0B_off.mat'};
tic 
q = animalBehaviorAnalysis('PROJECT_NAME',projectName,...
    'SAVE_LOCATION',savepath,...
    'MERGE', mergelist);
q.computeAnimalWiseProperties;
toc
saveFile = fullfile(q.meta.save_location,[q.meta.project_name '.mat']);
save(saveFile,'q');

%%

cd 'Z:\721430-MIMED_Geriatrics-DMIC06\Biostat\vanderwyk\projects\hoffman\data'


%%




%%
            figure;
            hold on;
            title(inputs.VAR,'Interpreter', 'none');
            boxplot(obj.statistics.uni.(inputs.VAR),obj.statistics.uni.group_label2);
            saveas(gcf,fullfile(obj.meta.save_location, 'oneway',inputs.TITLE, [inputs.VAR '.fig']));
            hold off;
            close(gcf);






