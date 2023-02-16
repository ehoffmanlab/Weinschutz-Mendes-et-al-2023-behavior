classdef animalBehaviorAnalysis < handle


    properties
        meta;
        %         save_location;
        %         project_name;

        files;
        %         startup_file;
        %         grouping_file;
        %         viewpoint_file;
        %         timing_file;
        %         merged_file;

        data;
        %         G;
        %         viewpoint_data;
        %         timing_data;
        %         merged_data;
        processed ;
        results;
        statistics ;

    end

    properties (Access = protected)
        %         Data ;
        status;


    end

    methods
        function obj = animalBehaviorAnalysis(varargin)
            % ANIMALBEHAVIORANALYSIS - Constructor
            %
            % USAGE:
            % O = animalBehaviorAnalysis('PROJECT_NAME',P,...
            %                            'SAVE_LOCATION',S,...
            %                            'ACCLIMATION_TIMING', A,...
            %                            'EVENT_TIMING',ET,...
            %                            'BASELINE_TIMING',BT,...
            %                            'GROUPING_FILE', G,...
            %                            'VIEWPOINT_FILE',V,...
            %                            'TIMING_FILE',T);
            %
            % ------------------------------------------------------------
            % REQUIRED ARGUMENTS:
            % 'PROJECT_NAME',P -- string, project name
            % 'SAVE_LOCATION',S -- string, root save location.  Files will
            %        be saved under fullfile(SAVE_LOCATION,PROJECT_NAME)
            % 'ACCLIMATION_TIMING', A -- 1x2 numeric array delineating the
            %        start and end of acclimation period in seconds
            % 'EVENT_TIMING',ET -- 1x3 numeric array delineating the
            %        pre-event period, the event period, and the post-event period in seconds
            % 'BASELINE_TIMING',BT -- 1x2 numeric array delineating the
            %        start and end of baseline period in seconds
            % 'GROUPING_FILE', G  -- cell array of full path to grouping
            %           files
            % 'VIEWPOINT_FILE',V -- cell array of full path to viewpoint
            %           files
            % 'TIMING_FILE',T -- string of full path to SINGLE timing file
            %
            %
            % OPTIONAL ARGUMENTS:
            % OUTPUTS:
            % O -- and animalBehaviorAnalysis object
            % NOTES:
            %  1) Currently very little by way of file verification
            %  2) Only a single common experimental timing is permitted
            %  3) Although a pre-event period is requested, no operations
            %  are currently implemented in it.
            %  4) See file requirements file for details on the assumed
            % format for timing, grouping, and viewpoint files.


            obj.meta.version = '2.0' ;
            obj.initArt(obj.meta.version);

            obj.status = 'INTIALIZING...';

            obj.meta = struct('save_location',[],...
                'project_name',[]);
            obj.files = struct('startup_file',[],...
                'grouping_file',[],...
                'viewpoint_file',[],...
                'timing_file',[],...
                'merged_file',[]);
            obj.data = struct('G',[],...
                'viewpoint_data',[],...
                'timing_data',[],...
                'merged_data',[]);
            obj.results = struct();

            %intialize some variables
            validLogical = [0 1];


            % components need for graphing
            path(path,fullfile(fileparts(which(mfilename)),'BL','boundedline'))
            path(path,fullfile(fileparts(which(mfilename)),'BL','catuneven'))
            path(path,fullfile(fileparts(which(mfilename)),'BL','Inpaint_nans'))
            path(path,fullfile(fileparts(which(mfilename)),'BL','singlepatch'))


            obj.status = 'PARSING COMMAND LINE ARGS...';
            % now parse any command line arguments
            try
                clap = inputParser; % Command Line Argument Parser
                clap.KeepUnmatched = 1;
                clap.addParameter('PROJECT_NAME', ...
                    [],...
                    @(x) ischar(x));

                clap.addParameter('SAVE_LOCATION', ...
                    [],...
                    @(x) ischar(x));

                clap.addParameter('GROUPING_FILE', ...
                    [],...
                    @(x) iscellstr(x) || isempty(x) || ischar(x));
                clap.addParameter('VIEWPOINT_FILE', ...
                    [],...
                    @(x) iscellstr(x) || isempty(x) || ischar(x));
                clap.addParameter('TIMING_FILE', ...
                    [],...
                    @(x) iscellstr(x) || isempty(x) || ischar(x));
                clap.addParameter('RECORDING_RESOLUTION', ...
                    1,...
                    @(x) isnumeric(x));
             
                clap.addParameter('ACCLIMATION_TIMING', ...
                    [],...
                    @(x) isnumeric(x));
                clap.addParameter('EVENT_TIMING', ...
                    [],...
                    @(x) isnumeric(x));
                
                clap.addParameter('BASELINE_TIMING', ...
                    [],...
                    @(x) isnumeric(x));
                clap.addParameter('EVENT_INDEX', ...
                    [],...
                    @(x) ischar(x));
                clap.addParameter('EVENT_LABEL', ...
                    [],...
                    @(x) ischar(x) || isempty(x) );
                clap.addParameter('DEP_VAR', ...
                    [],...
                    @(x) ischar(x));
                clap.addParameter('OVERWRITE_ON_LOAD', ...
                    0,...
                    @(x) isnumeric(x));
                clap.addParameter('MERGE', ...
                    {},...
                    @(x) iscell(x));
                clap.parse(varargin{:});
                clap_inputs = clap.Results;
            catch err
                obj.status = sprintf('FUNCTION "%s" CALL ERROR: %s', upper(obj.current_function),err.message);
                obj.help('animalBehaviorAnalysis');
                return;
            end

            %Check command line arguments

            %query user if unspecified
            if isempty(clap_inputs.PROJECT_NAME)
                clap_inputs.PROJECT_NAME = input('\nEnter Project Name:     ','s');
            end

            if isempty(clap_inputs.SAVE_LOCATION)
                clap_inputs.SAVE_LOCATION = uigetdir(pwd,'Select Analysis Save Location');

            end
            obj.meta.rootAnalysisOutput = clap_inputs.SAVE_LOCATION;
            %

            % Load and save timing info
            if ~isempty(clap_inputs.ACCLIMATION_TIMING)
                obj.meta.acclimationPeriodStart = clap_inputs.ACCLIMATION_TIMING(1);
                obj.meta.acclimationPeriodEnd= clap_inputs.ACCLIMATION_TIMING(2);
                obj.meta.acclimation.pre_event_period = 0;
                obj.meta.acclimation.event_period = 1;
                obj.meta.acclimation.post_event_period  = 0;
            end
            if ~isempty(clap_inputs.BASELINE_TIMING)
                obj.meta.commonBaselinePeriodStart = clap_inputs.BASELINE_TIMING(1);
                obj.meta.commonBaselinePeriodEnd= clap_inputs.BASELINE_TIMING(2);
                obj.meta.common_baseline.pre_event_period = 0;
                obj.meta.common_baseline.event_period = 1;
                obj.meta.common_baseline.post_event_period  = 0;
            end
            obj.meta.event = struct();
            if ~isempty(clap_inputs.EVENT_TIMING)
                obj.meta.event.pre_event_period = clap_inputs.EVENT_TIMING(1);
                obj.meta.event.event_period = clap_inputs.EVENT_TIMING(2);
                obj.meta.event.post_event_period  = clap_inputs.EVENT_TIMING(3);
            else
                obj.meta.event.pre_event_period = [];
                obj.meta.event.event_period = [];
                obj.meta.event.post_event_period  = [];

            end

            obj.meta.recordingResolution = clap_inputs.RECORDING_RESOLUTION ;
            
            if ~isempty(clap_inputs.EVENT_INDEX)
                obj.meta.event.event_index = clap_inputs.EVENT_INDEX;
            else
                obj.meta.event.event_index = [];
            end
            if ~isempty(clap_inputs.EVENT_LABEL)
                obj.meta.event.event_label = clap_inputs.EVENT_LABEL;
            else
                obj.meta.event.event_label = [];
            end
            if ~isempty(clap_inputs.DEP_VAR)
                obj.meta.dep_var = clap_inputs.DEP_VAR;
            else
                obj.meta.dep_var = [];
            end


            %             % Load commands if not specified -- this has not been tested in
            %             % a while.
            %             if isempty(clap_inputs.GROUPING_FILE)
            %                 clap_inputs.GROUPING_FILE = uigetdir(pwd,'Select GROUPING file');
            %             end
            %
            %             if isempty(clap_inputs.VIEWPOINT_FILE)
            %                 clap_inputs.VIEWPOINT_FILE = uigetdir(pwd,'Select VIEWPOINT file');
            %             end
            %
            %
            %             if isempty(clap_inputs.TIMING_FILE)
            %                 clap_inputs.TIMING_FILE = uigetdir(pwd,'Select TIMING file');
            %             end

            if isempty(clap_inputs.GROUPING_FILE)
                obj.files.grouping_file = [] ;
            else
                obj.files.grouping_file = clap_inputs.GROUPING_FILE;
            end
            if isempty(clap_inputs.VIEWPOINT_FILE)
                obj.files.viewpoint_file = [] ;
            else
                obj.files.viewpoint_file = clap_inputs.VIEWPOINT_FILE;
            end
            if isempty(clap_inputs.TIMING_FILE)
                obj.files.timing_file = [] ;
            else
                obj.files.timing_file = clap_inputs.TIMING_FILE;
            end

            %make sure save location exists
            obj.meta.save_location = fullfile(clap_inputs.SAVE_LOCATION,clap_inputs.PROJECT_NAME);
            if exist(fullfile(obj.meta.rootAnalysisOutput),'dir')
                if ~exist(obj.meta.save_location,'dir')
                    try
                        mkdir(obj.meta.save_location);
                    catch err
                        obj.status = sprintf('FUNCTION "%s" CALL ERROR: %s', upper(obj.current_function),err.message);
                        return;
                    end
                end
            else
                obj.status = 'ERROR: Specified save location does not exist';
                return;
            end
            obj.status = 'COMMAND LINE ARGS PARSED!';

            %             obj.meta.version = 2.0 ;


            %assign inputs
            %             obj.files.startup_file = sap_inputs.STARTUP_FILE;

            obj.meta.project_name = clap_inputs.PROJECT_NAME;



            obj.files.result_file = fullfile(obj.meta.rootAnalysisOutput,...
                obj.meta.project_name,...
                'result.txt');

            %             obj.meta.input_location = fileparts(obj.files.viewpoint_file{1}) ;


            %             saved_merge = exist(fullfile(obj.meta.save_location,'merged_data_table.mat'),'file');

            %             if saved_merge == 2 && clap_inputs.OVERWRITE_ON_LOAD == 0
            %                 obj.status = 'LOADING PRIOR MERGED DATA...';
            %                 temp = load(fullfile(obj.meta.save_location,'merged_data_table.mat'));
            %                 obj.data.merged_data = temp.t ;
            %                 obj.halt_xls_load = 1 ;
            %                 obj.status = 'DONE!';
            %             else
            %                 obj.halt_xls_load = 0 ;
            %             end
            % Load files


            if isempty(clap_inputs.MERGE)
                obj.status = 'LOADING FILES...';

                s = obj.loadFiles('FILE_SET','GROUPING');
                if s.error
                    obj.meta.loadError = 'GROUPING LOAD ERROR';
                    return;
                end
                s = obj.loadFiles('FILE_SET','VIEWPOINT');
                if s.error
                    obj.meta.loadError = 'VIEWPOINT LOAD ERROR';
                    return;
                end
                s = obj.loadFiles('FILE_SET','TIMING');
                if s.error
                    obj.meta.loadError = 'TIMING LOAD ERROR';
                    return;
                end

                obj.meta.loadError = [];

                obj.status = 'FILES LOADED!';

                % Merge input files
                m = obj.mergeDataFiles();

                if m.error
                    return;
                end

                % Do a basic validity check on the data before proceeding 
                obj.validityCheck();
                
                
                %split data into numeric and non-numeric types
                obj.splitData ;
                % overide preperiod if 0
                if obj.meta.event.pre_event_period == 0
                    obj.meta.event.pre_event_period = 1 ;

                end


                obj.status = 'LABEL EVENTS!';

                obj.labelEvents ;

                obj.status = 'DONE!';

                %             t = obj.data.merged_data;
                %             save(fullfile( obj.meta.save_location,'merged_data_table.mat'),'t');


                obj.status = 'RESTRUCTURE DATA!';

                obj.restructureData ;

                obj.status = 'DONE!';
                obj.status = 'PROCESS INDIVIDUAL DATA!';

                obj.processIndividualData ;

                obj.status = 'DONE!';

            else
                obj.status = 'MERGING FILES...';
                obj.mergeobjs('PROJECT_NAME',clap_inputs.PROJECT_NAME,...
                    'SAVE_LOCATION',clap_inputs.SAVE_LOCATION,...
                    'MERGE', clap_inputs.MERGE);
                obj.status = 'DONE!';

            end


            obj.status = 'PROCESS GROUPWISE DATA!';

            obj.processGroupwiseData ;

            obj.status = 'DONE!';

            %             str = input(prompt,'s')
            %             folder_name = uigetdir(start_path)
        end
        function set.status(obj,varargin)
            % FUNCTION SET.STATUS (s)
            %
            % USAGE:
            % obj.set.()
            %
            % ------------------------------------------------------------
            % REQUIRED ARGUMENTS:
            %
            % OPTIONAL ARGUMENTS:
            %
            % OUTPUTS:
            %
            % NOTES:
            %


            %WIERD BEHAVIOR ON LOAD
            try
                %check inputs
                p = inputParser;
                p.addRequired('value', ...
                    @(x) ischar(x));
                p.parse(varargin{:});
                inputs = p.Results;
            catch err
                %fprintf('\nSTATUS ERROR: %s\n\n', err.message);
                return;
            end
            obj.status = upper(inputs.value);
            fprintf('\nSTATUS: %s\n',obj.status);
            obj.status = [];
        end
        function mergeobjs(obj,varargin)
            try
                clap = inputParser; % Command Line Argument Parser
                clap.KeepUnmatched = 1;
                clap.addParameter('PROJECT_NAME', ...
                    [],...
                    @(x) ischar(x));

                clap.addParameter('SAVE_LOCATION', ...
                    [],...
                    @(x) ischar(x));
                clap.addParameter('MERGE', ...
                    {},...
                    @(x) iscell(x));
                clap.parse(varargin{:});
                clap_inputs = clap.Results;
            catch err
                obj.status = sprintf('FUNCTION "%s" CALL ERROR: %s', upper(obj.current_function),err.message);
                obj.help('merge');
                return;
            end

            % check version and ABA class
            % not implemented yet

            % absorb the first obj
            x = load(clap_inputs.MERGE{1});
            obj.meta = x.n.meta ;
            obj.files = x.n.files ;
            obj.data = x.n.data ;
            obj.processed = x.n.processed ;

            for o = 2 : numel(clap_inputs.MERGE)


                %max animal_id in current data
                ID_inc = max(obj.data.group_data.animal_id) ;


                x = load(clap_inputs.MERGE{o});
                %reconcile files

                % check timing of potential loaded file 
                curr_timing = size(obj.processed.block{1}.ind_timeseries,2) ;
                merge_timing = size(x.n.processed.block{1}.ind_timeseries,2);
                
                if curr_timing ~= merge_timing 
                    errorStruct.message = sprintf('During OBJ merge %s, timing mismatch was detected in %s',obj.meta.project_name, x.n.meta.project_name ) ;
                    errorStruct.identifier = 'mergeobjs:mergeError';
                    error(errorStruct)
                end

                %reconcile meta.group_cross
                newGroup = ['groupNum_', sprintf('%d',o) ] ;
                x.n.meta.group_cross.Properties.VariableNames = {newGroup,'group_label2'};


                %                 nr.Properties.VariableNames = {'groupNum_2', 'group_label2'};
                %                 x.n.meta.group_cross = [x.n.meta.group_cross;nr];
                %
                C = outerjoin(obj.meta.group_cross,x.n.meta.group_cross,'Keys','group_label2') ;
                C.group_label2 = repmat({''},height(C),1);

                newGroupNum = max(C.groupNum) + 1 ;
                for r = 1 :height(C)
                    if isnan(C.groupNum(r))
                        C.groupNum(r) = newGroupNum ;
                        newGroupNum = newGroupNum + 1 ;

                    end

                    if isempty(C.group_label2_left{r})
                        C.group_label2{r} = C.group_label2_right{r} ;
                    else
                        C.group_label2{r} = C.group_label2_left{r} ;
                    end


                end

                C.group_label2_left = [] ;
                C.group_label2_right = [] ;
                sortrows(C,'groupNum') ;
                obj.meta.group_cross = C ;


                %data to be reconciled

                x.n.data.group_data.animal_id = x.n.data.group_data.animal_id + ID_inc ;
                x.n.data.group_data.groupNum;
                for r = 1 : size(x.n.data.group_data.groupNum,1)
                    x.n.data.group_data.groupNum(r) = obj.meta.group_cross.groupNum(x.n.data.group_data.groupNum(r)==obj.meta.group_cross.(newGroup));
                end
                obj.data.group_data = vertcat(obj.data.group_data,x.n.data.group_data) ;


                %per block
                for b = 1 : numel(obj.processed.block)
                    x.n.processed.block{b}.ind_desc(:,x.n.meta.PDVariableNames.animal_id) = x.n.processed.block{b}.ind_desc(:,x.n.meta.PDVariableNames.animal_id)  + ID_inc ;
                    try
                        for r = 1 : size(x.n.processed.block{b}.ind_desc(:,x.n.meta.PDVariableNames.groupNum),1)
                            x.n.processed.block{b}.ind_desc(r,x.n.meta.PDVariableNames.groupNum) = obj.meta.group_cross.groupNum(x.n.processed.block{b}.ind_desc(r,x.n.meta.PDVariableNames.groupNum)==obj.meta.group_cross.(newGroup));
                        end
                    catch err
                    end
                    obj.processed.block{b}.ind_desc =    vertcat(obj.processed.block{b}.ind_desc, x.n.processed.block{b}.ind_desc) ;


                    % Currently the animal id column is hardcoded for the
                    % timeseries data
                    try
                        x.n.processed.block{b}.ind_timeseries(:,1) = x.n.processed.block{b}.ind_timeseries(:,1) + ID_inc;
                        for r = 1 :size( x.n.processed.block{b}.ind_timeseries,1)
                            x.n.processed.block{b}.ind_timeseries(r,2) = obj.meta.group_cross.groupNum(x.n.processed.block{b}.ind_timeseries(r,2)==obj.meta.group_cross.(newGroup));
                        end
                    catch err
                    end
                    try 
                    obj.processed.block{b}.ind_timeseries =    vertcat(obj.processed.block{b}.ind_timeseries, x.n.processed.block{b}.ind_timeseries) ;
                    catch err 
                        pause ;
                    end 
                end

                % clear group descriptives ... this will need to be
                % recomputed

                obj.processed.block{b}.grp_desc = [];




            end

            obj.meta.save_location = fullfile(clap_inputs.SAVE_LOCATION,clap_inputs.PROJECT_NAME);
            obj.meta.project_name = clap_inputs.PROJECT_NAME ;

            obj.results = [];

        end
        function validityCheck(obj)
            
            
            % check time resolution
             if (obj.data.viewpoint_data.xEnd(1) - obj.data.viewpoint_data.start(1)) ~= obj.meta.recordingResolution
                 errorStruct.message = sprintf('Recording time resolution from viewpoint (%d) does not match command (%d).',(obj.data.viewpoint_data.xEnd(1) - obj.data.viewpoint_data.start(1)),obj.meta.recordingResolution) ;
                 errorStruct.identifier = 'validityCheck:viewpoint';
                 error(errorStruct)
             end
            
        end
        function processIndividualData(obj, varargin)
            dv = obj.meta.NVariableNames.(obj.meta.dep_var);
            animal_id = obj.meta.NVariableNames.animal_id ;
            groupNum = obj.meta.NVariableNames.groupNum ;
            %             build stat meta struct
            obj.meta.PDVariableNames = struct('count', 0);
            obj.registerNewDataColumn('PDVariableNames', 'animal_id');
            obj.registerNewDataColumn('PDVariableNames', 'groupNum');
            obj.registerNewDataColumn('PDVariableNames', 'mean');
            obj.registerNewDataColumn('PDVariableNames', 'std');
            obj.registerNewDataColumn('PDVariableNames', 'conmean');
            obj.registerNewDataColumn('PDVariableNames', 'propmv');
            obj.registerNewDataColumn('PDVariableNames', 'propNum');




            obj.processed.block = [];

            nBlocks = numel(obj.data.nbs);
            for b = 1 : nBlocks
                %                 thisBlock = obj.data.nbs{b}.stacked ;


                try

                thisMean = nanmean(obj.data.nbs{b}.stacked(:,:,dv),2);
                thisStd = nanstd(obj.data.nbs{b}.stacked(:,:,dv),0,2);
                tempData = obj.data.nbs{b}.stacked(:,:,dv) ;
                zeroEntries = find(tempData(:,:)==0) ;
                tempData(zeroEntries) = NaN ;
                thisConMean = nanmean(tempData(:,:),2) ;
                thisPropNum = sum(obj.data.nbs{b}.stacked(:,:,dv)>0,2) ;
                thisPropDenom = size(obj.data.nbs{b}.stacked(:,:,dv),2) ;
                thisPropMv = thisPropNum / thisPropDenom ;
                thisAID = obj.data.nbs{b}.stacked(:,1,animal_id) ;
                thisGID = obj.data.nbs{b}.stacked(:,1,groupNum) ;

                thisTime = obj.data.nbs{b}.stacked(:,:,dv);
                catch err
                    
                end
                
                obj.processed.block{b} = struct('name', obj.data.nbs{b}.name,...
                    'ind_desc', [thisAID,...
                    thisGID,...
                    thisMean,...
                    thisStd,...
                    thisConMean,...
                    thisPropMv,...
                    thisPropNum],...
                    'ind_timeseries',[thisAID,...
                    thisGID,...
                    thisTime]);



            end
            obj.computeDerivedVariables() ;




        end
        function processGroupwiseData(obj,varargin)
            groupNum = obj.meta.PDVariableNames.groupNum;
            whichstats = {'mean','sem','std','min','max','meanci','numel'};

            nBlocks = numel(obj.processed.block) ;
            for n = 1 : nBlocks
                [a,b,c,d,e,f,g] = grpstats(obj.processed.block{n}.ind_desc,obj.processed.block{n}.ind_desc(:,groupNum),whichstats);
                obj.processed.block{n}.grp_desc.mean = a ;
                obj.processed.block{n}.grp_desc.sem = b ;
                obj.processed.block{n}.grp_desc.std = c ;
                obj.processed.block{n}.grp_desc.min = d ;
                obj.processed.block{n}.grp_desc.max = e ;
                obj.processed.block{n}.grp_desc.lci = f(:,:,1) ;
                obj.processed.block{n}.grp_desc.uci = f(:,:,2) ;
                obj.processed.block{n}.grp_desc.numel = f(:,:,2) ;


            end
        end

        function stats = computeAnimalWiseProperties(obj, varargin)
            stats = [] ;

            groupNum = obj.meta.PDVariableNames.groupNum;
            animal_id = obj.meta.PDVariableNames.animal_id;
            c_mean = obj.meta.PDVariableNames.mean;
            c_std = obj.meta.PDVariableNames.std;
            c_conmean = obj.meta.PDVariableNames.conmean;
            c_propmv = obj.meta.PDVariableNames.propmv;
            c_propNum = obj.meta.PDVariableNames.propNum;
            c_nrmMean = obj.meta.PDVariableNames.nrmMean;
            c_nrmConMean = obj.meta.PDVariableNames.nrmConMean;
            c_nrmPropMv = obj.meta.PDVariableNames.nrmPropMv;

            %set up the results
            %             obj.meta.ResultVariableNames = struct('count', 0);
            %             obj.registerNewDataColumn('ResultVariableNames', 'animal_id');
            %             obj.registerNewDataColumn('ResultVariableNames', 'groupNum');
            root_res_dir = fullfile(obj.meta.save_location,'results');
            %             if ~exist(root_res_dir,7)
            %                    mkdir(root_res_dir) ;
            %             end

            obj.results.block = [] ;

            if ~exist(fullfile(obj.meta.save_location,'results'),'dir')
               mkdir(fullfile(obj.meta.save_location,'results')) ;
            end
             if ~exist(fullfile(obj.meta.save_location,'results','oneway'),'dir')
               mkdir(fullfile(obj.meta.save_location,'results','oneway')) ;
            end
            for b = 1 : numel(obj.processed.block)
                %create result directory

                if ~exist(fullfile(obj.meta.save_location,'results','oneway',obj.processed.block{b}.name),'dir')
                    mkdir(fullfile(obj.meta.save_location,'results','oneway',obj.processed.block{b}.name)) ;
                end

                res_dir_block = fullfile(root_res_dir,obj.processed.block{b}.name);

                obj.results.block{b}.name = obj.processed.block{b}.name ;
                %%%%% stats for means %%%%%%%

                obj.results.block{b}.var_mean = struct()     ;
                obj.results.block{b}.var_mean.oneway = obj.onewayANOVA('DV',obj.processed.block{b}.ind_desc(:,c_mean) ,...
                    'GRP',obj.processed.block{b}.ind_desc(:,groupNum) );
%                 obj.zscore('DV',obj.processed.block{b}.ind_desc(:,c_mean) ,...
%                     'GRP',obj.processed.block{b}.ind_desc(:,groupNum),...
%                     'PROP', 0)
%                 % stats for std
                obj.results.block{b}.var_std = struct()     ;
                obj.results.block{b}.var_std.oneway = obj.onewayANOVA('DV',obj.processed.block{b}.ind_desc(:,c_std) ,...
                    'GRP',obj.processed.block{b}.ind_desc(:,groupNum) );

                % stats for conmean
                obj.results.block{b}.var_conmean = struct()     ;
                obj.results.block{b}.var_conmean.oneway = obj.onewayANOVA('DV',obj.processed.block{b}.ind_desc(:,c_conmean) ,...
                    'GRP',obj.processed.block{b}.ind_desc(:,groupNum) );

                % stats for propmv
                obj.results.block{b}.var_propmv = struct()     ;
                obj.results.block{b}.var_propmv.oneway = obj.onewayANOVA('DV',obj.processed.block{b}.ind_desc(:,c_propmv) ,...
                    'GRP',obj.processed.block{b}.ind_desc(:,groupNum) );

                % stats for nrmMean
                obj.results.block{b}.var_nrmMean = struct()     ;
                obj.results.block{b}.var_nrmMean.oneway = obj.onewayANOVA('DV',obj.processed.block{b}.ind_desc(:,c_nrmMean) ,...
                    'GRP',obj.processed.block{b}.ind_desc(:,groupNum) );

                % stats for nrmConMean
                obj.results.block{b}.var_nrmConMean = struct()     ;
                obj.results.block{b}.var_nrmConMean.oneway = obj.onewayANOVA('DV',obj.processed.block{b}.ind_desc(:,c_nrmConMean) ,...
                    'GRP',obj.processed.block{b}.ind_desc(:,groupNum) );

                % stats for nrmPropMv
                obj.results.block{b}.var_nrmPropMv = struct()     ;
                obj.results.block{b}.var_nrmPropMv.oneway = obj.onewayANOVA('DV',obj.processed.block{b}.ind_desc(:,c_nrmPropMv) ,...
                    'GRP',obj.processed.block{b}.ind_desc(:,groupNum) );


            end
            obj.printANOVAtables ;
            obj.printGroupTable ;
            obj.simpleBoxPlots ;
            obj.printPsummaryTable() ;
            obj.printIndividualTable() ;
            obj.printTimeSeriesGraphs() ;

        end
        function [l_data, w_data] = pullData(obj,varargin) 
            
            %this might be able to be streamlined a bit
            
            l_data = [];
            w_data = [];
            arg = inputParser; % Command Line Argument Parser
            
            arg.addParameter('SCOPE', ...
                'GROUP',...
                @(x) ischar(x));
            
            arg.parse(varargin{:});
            arg_inputs = arg.Results;
            
           
            if strcmpi('GROUP',arg_inputs.SCOPE)
                for b = 1 : numel(obj.processed.block)
                    
                    t = array2table(obj.processed.block{b}.grp_desc.mean);
                    t.Properties.VariableNames =  {'animal_id' 'groupNum' 'mean' 'std' 'conmean' 'propmv' 'propNum' 'nrmMean' 'nrmConMean' 'nrmPropMv'};
                    t.animal_id = [];
                    t.groupLabel = obj.returnGroups(t.groupNum) ;
                    block_name =  obj.results.block{b}.name ;
                    t.block =   repmat({block_name},height(t),1) ;
                    
                    %reorder columns
                    oldvariables = t.Properties.VariableNames;
                    newvariables = {'groupLabel' 'groupNum' 'block' 'mean' 'std' 'conmean' 'propmv' 'propNum' 'nrmMean' 'nrmConMean' 'nrmPropMv'};
                    [~,LOCB] = ismember(newvariables,oldvariables);
                    t = t(:,LOCB);

                    
                    wt = t ;
                    wt.block = [];
                    
                    

                    
                    for v = 3 : numel(wt.Properties.VariableNames)
                        wt.Properties.VariableDescriptions{v} = [block_name '(' wt.Properties.VariableNames{v}  ')']; 
                        wt.Properties.VariableNames{v} =  strrep(strrep(strrep(wt.Properties.VariableDescriptions{v},'(','_'),')',''),'-','m');
                    end
                    
                    
                    if b == 1
                        l_data = t ;
                        w_data = wt ;
                    else
                        l_data = [l_data ; t] ;
                        wt.groupNum = [];
                        wt.groupLabel = [];
                        w_data = [w_data wt] ;
                    end
                end
                
            else
                for b = 1 : numel(obj.processed.block)
                    
                    t = array2table(obj.processed.block{b}.ind_desc);
                    t.Properties.VariableNames =  {'animal_id' 'groupNum' 'mean' 'std' 'conmean' 'propmv' 'propNum' 'nrmMean' 'nrmConMean' 'nrmPropMv'};
                   
                    t.groupLabel = obj.returnGroups(t.groupNum) ;
                    
                    block_name =  obj.results.block{b}.name ;
                    t.block =   repmat({block_name},height(t),1) ;
                    
                    
                     %reorder columns
                    oldvariables = t.Properties.VariableNames;
                    newvariables = {'animal_id' 'groupLabel' 'groupNum' 'block' 'mean' 'std' 'conmean' 'propmv' 'propNum' 'nrmMean' 'nrmConMean' 'nrmPropMv'};
                    [~,LOCB] = ismember(newvariables,oldvariables);
                    t = t(:,LOCB);
                    
                     wt = t ;
                    wt.block = [];
                    
                     for v = 4 : numel(wt.Properties.VariableNames)
                        wt.Properties.VariableDescriptions{v} = [block_name '(' wt.Properties.VariableNames{v}  ')']; 
                        wt.Properties.VariableNames{v} =  strrep(strrep(strrep(wt.Properties.VariableDescriptions{v},'(','_'),')',''),'-','m');
                     end
                    
                    if b == 1
                        l_data = t ;
                        w_data = wt ;
                    else
                        l_data = [l_data ; t] ;
                        wt.groupNum = [];
                        wt.groupLabel = [];
                        wt.animal_id = [];
                        w_data = [w_data wt] ;
                    end
                end
            end
            
            
            
            
        end
        function printANOVAtables(obj)
            animal_id = obj.meta.PDVariableNames.animal_id;
            groupNum = obj.meta.PDVariableNames.groupNum;
            c_mean = obj.meta.PDVariableNames.mean;
            c_std = obj.meta.PDVariableNames.std;
            c_conmean = obj.meta.PDVariableNames.conmean;
            c_propmv = obj.meta.PDVariableNames.propmv;
            c_propNum = obj.meta.PDVariableNames.propNum;
            c_nrmMean = obj.meta.PDVariableNames.nrmMean;
            c_nrmConMean = obj.meta.PDVariableNames.nrmConMean;
            c_nrmPropMv = obj.meta.PDVariableNames.nrmPropMv;

            
            for b = 1 : numel(obj.results.block)
                rtab = cell2table(obj.results.block{b}.var_mean.oneway.tbl(2:end,2:end));
                rtab.Properties.VariableNames = genvarname(obj.results.block{b}.var_mean.oneway.tbl(1,2:end)) ;
                rtab.Properties.VariableNames{end} = 'P' ;

                rtab.Properties.RowNames = genvarname(obj.results.block{b}.var_mean.oneway.tbl(2:end,1)');
                    delete(fullfile( obj.meta.save_location, 'results','oneway',obj.processed.block{b}.name,'anova_table.csv' ))
                    delete(fullfile( obj.meta.save_location, 'results','oneway',obj.processed.block{b}.name,'anova_table.xlsx' ))
                    writetable(rtab, fullfile( obj.meta.save_location, 'results','oneway',obj.processed.block{b}.name,'anova_table.csv' ) ,'WriteRowNames' , 1);

                if ~isempty(obj.results.block{b}.var_mean.oneway.ctab)
                   delete(fullfile( obj.meta.save_location, 'results','oneway',obj.processed.block{b}.name,'multComp_table.xlsx' ))
                   delete(fullfile( obj.meta.save_location, 'results','oneway',obj.processed.block{b}.name,'multComp_table.csv' ))

                    writetable(obj.results.block{b}.var_mean.oneway.ctab, fullfile( obj.meta.save_location, 'results','oneway',obj.processed.block{b}.name,'multComp_table.csv' ) ,'WriteRowNames' , 1);
                end

            end
        end
        function printPsummaryTable(obj)
            animal_id = obj.meta.PDVariableNames.animal_id;
            groupNum = obj.meta.PDVariableNames.groupNum;
            c_mean = obj.meta.PDVariableNames.mean;
            c_std = obj.meta.PDVariableNames.std;
            c_conmean = obj.meta.PDVariableNames.conmean;
            c_propmv = obj.meta.PDVariableNames.propmv;
            c_propNum = obj.meta.PDVariableNames.propNum;
            c_nrmMean = obj.meta.PDVariableNames.nrmMean;
            c_nrmConMean = obj.meta.PDVariableNames.nrmConMean;
            c_nrmPropMv = obj.meta.PDVariableNames.nrmPropMv;

            for b = 1 : numel(obj.results.block)
                  vars = {'mean' 'std' 'conmean' 'propmv' 'nrmMean' 'nrmConMean' 'nrmPropMv'} ;
                  pvals = nan(numel(vars),1);
                  for v = 1:numel(vars)
                      try
                          pvals(v) = obj.results.block{b}.(['var_' vars{v}]).oneway.p;
                      catch
                      end
                  end
                  block_name =  repmat({obj.results.block{b}.name},numel(vars),1) ;
                  descr = repmat({'testing significance of group'},numel(vars),1) ;
                  t = table(block_name,vars',pvals,descr ) ;
                  t.Properties.VariableNames =  {'Block' 'Variable' 'PValue' 'Note'};
%                   writetable(t, fullfile( obj.meta.save_location, 'results', obj.processed.block{b}.name,'pvalue_summary.xlsx')  );

                  if b == 1
                      agg_table = t ;
                  else
                      agg_table = [agg_table ; t] ;
                  end

            end
            writetable(agg_table, fullfile( obj.meta.save_location, 'results', 'pvalue_summary_mean.xlsx')  );

        end
        function printTimeSeriesGraphs(obj)

%             Full time series un-normed
            groups = obj.meta.group_cross.groupNum ;
            figure;
            hold on;
            title('Full Raw Timeseries','Interpreter', 'none');
            cmap = lines(numel(groups));
            for g = numel(groups):-1:1
                %compute timeise SEM per group
                gn = groups(g) ;

                q = obj.processed.block{1}.ind_timeseries(obj.processed.block{1}.ind_timeseries(:,2)==gn,:);
                y = mean(q(:,3:end));
                x = 1:numel(y);
                sem = std(q(:,3:end))/sqrt(size(q(:,3:end),1)) ;
                boundedline(x,y,sem,'alpha','cmap', cmap(g,:), 'transparency', 0.5);
            end
            h = get(gca,'Children');
            h1 = h(1:2:end);
            legend(h1,obj.returnGroups(groups),'Interpreter','none');

            saveas(gcf,fullfile(obj.meta.save_location,'results','Full_Raw_Timeseries.fig'));
            hold off;
            close(gcf);


%             No baseline raw timeseries
             groups = obj.meta.group_cross.groupNum ;
            figure;
            hold on;
            title('No Baseline Raw Timeseries','Interpreter', 'none');
            cmap = lines(numel(groups));
%           we are picking a starting point for the "event" time periods
%           that is just before the start of the actual events
            event_starts = obj.meta.commonBaselinePeriodEnd - obj.meta.commonBaselinePeriodStart - 5 ;
            for g = numel(groups):-1:1
                %compute timeise SEM per group
                gn = groups(g) ;

                q = obj.processed.block{1}.ind_timeseries(obj.processed.block{1}.ind_timeseries(:,2)==gn,:);
                y = mean(q(:,event_starts:end));
                x = event_starts:numel(y)+event_starts-1;
                sem = std(q(:,event_starts:end))/sqrt(size(q(:,event_starts:end),1)) ;
                boundedline(x,y,sem,'alpha','cmap', cmap(g,:), 'transparency', 0.5);
            end
            h = get(gca,'Children');
            h1 = h(1:2:end);
            legend(h1,obj.returnGroups(groups),'Interpreter','none');

            saveas(gcf,fullfile(obj.meta.save_location,'results','NoBaseline_Raw_Timeseries.fig'));
            hold off;
            close(gcf);


            % No baseline, proportion move
            groups = obj.meta.group_cross.groupNum ;
            figure;
            hold on;
            title('No Baseline Proportion Moving Timeseries','Interpreter', 'none');
            cmap = lines(numel(groups));
%           we are picking a starting point for the "event" time periods
%           that is just before the start of the actual events
            event_starts = obj.meta.commonBaselinePeriodEnd - obj.meta.commonBaselinePeriodStart - 5 ;
            for g = numel(groups):-1:1
                %compute timeise SEM per group
                gn = groups(g) ;

                % we have id the proportion moving


                q = obj.processed.block{1}.ind_timeseries(obj.processed.block{1}.ind_timeseries(:,2)==gn,:);
                y = mean(q(:,event_starts:end)>0);
                x = event_starts:numel(y)+event_starts-1;
                sem = zeros(1,numel(x));
                plot(x,y);
            end
            h = get(gca,'Children');
            h1 = h(1:1:end);
            legend(h1,obj.returnGroups(groups),'Interpreter','none');

            saveas(gcf,fullfile(obj.meta.save_location,'results','NoBaseline_propmv_Timeseries.fig'));
            hold off;
            close(gcf);



        end
        function printIndividualTable(obj)
            
            %pull data 
            [l,w] = obj.pullData('SCOPE','IND');
            
                        
            %print wide format
            writetable(w, fullfile( obj.meta.save_location, 'results', 'individual_descriptives_wide.csv')  );
            
            %print long format
            writetable(l, fullfile( obj.meta.save_location, 'results', 'individual_descriptives_long.csv')  );
            
        end
        function printIndividualTable_arch(obj)
            animal_id = obj.meta.PDVariableNames.animal_id;
            groupNum = obj.meta.PDVariableNames.groupNum;
            c_mean = obj.meta.PDVariableNames.mean;
            c_std = obj.meta.PDVariableNames.std;
            c_conmean = obj.meta.PDVariableNames.conmean;
            c_propmv = obj.meta.PDVariableNames.propmv;
            c_propNum = obj.meta.PDVariableNames.propNum;
            c_nrmMean = obj.meta.PDVariableNames.nrmMean;
            c_nrmConMean = obj.meta.PDVariableNames.nrmConMean;
            c_nrmPropMv = obj.meta.PDVariableNames.nrmPropMv;

            for b = 1 : numel(obj.processed.block)

                  t = array2table(obj.processed.block{b}.ind_desc);
                  t.Properties.VariableNames =  {'animal_id' 'groupNum' 'mean' 'std' 'conmean' 'propmv' 'propNum' 'nrmMean' 'nrmConMean' 'nrmPropMv'};
                  t.animal_id = [];
                  t.groupLabel = obj.returnGroups(t.groupNum) ;

                  block_name =  repmat({obj.results.block{b}.name},height(t),1) ;
                  t.block = block_name ;

                  if b == 1
                      agg_table = t ;
                  else
                      agg_table = [agg_table ; t] ;
                  end
            end
            writetable(agg_table, fullfile( obj.meta.save_location, 'results', 'individual_descriptives.xlsx')  );

            % Now build and write a wide format 
            
            
        end
        function printGroupTable(obj)
            
            %pull data 
            [l,w] = obj.pullData('SCOPE','GROUP');
            
            %print wide format
            writetable(w, fullfile( obj.meta.save_location, 'results', 'group_descriptives_wide.csv')  );
            
            %print long format
            writetable(l, fullfile( obj.meta.save_location, 'results', 'group_descriptives_long.csv')  );
        end
        
        
        function printGroupTable_arch(obj)

            animal_id = obj.meta.PDVariableNames.animal_id;
            groupNum = obj.meta.PDVariableNames.groupNum;
            c_mean = obj.meta.PDVariableNames.mean;
            c_std = obj.meta.PDVariableNames.std;
            c_conmean = obj.meta.PDVariableNames.conmean;
            c_propmv = obj.meta.PDVariableNames.propmv;
            c_propNum = obj.meta.PDVariableNames.propNum;
            c_nrmMean = obj.meta.PDVariableNames.nrmMean;
            c_nrmConMean = obj.meta.PDVariableNames.nrmConMean;
            c_nrmPropMv = obj.meta.PDVariableNames.nrmPropMv;

            for b = 1 : numel(obj.processed.block)

                  t = array2table(obj.processed.block{b}.grp_desc.mean);
                  t.Properties.VariableNames =  {'animal_id' 'groupNum' 'mean' 'std' 'conmean' 'propmv' 'propNum' 'nrmMean' 'nrmConMean' 'nrmPropMv'};
                  t.animal_id = [];
                        t.groupLabel = obj.returnGroups(t.groupNum) ;

                  block_name =  obj.results.block{b}.name ;
                 
                  t.block =   repmat({block_name},height(t),1) ;

                  
                  wt = t ;
                  
                  wt.Properties.VariableDescriptions =  {'animal_id' ...
                      'groupNum' ...
                      [block_name '(mean)'] ...
                      [block_name '(std)'] ...
                      [block_name '(conmean)'] ...
                      [block_name '(propmv)'] ...
                      [block_name '(propNum)'] ...
                      [block_name '(nrmMean)'] ...
                      [block_name '(nrmConMean)'] ...
                      [block_name '(nrmPropMv)'] ... 
                      'block' };
                  for v = 1 : numel(wt.Properties.VariableNames)
                      wt.Properties.VariableNames{v} =  strrep(strrep(strrep(wt.Properties.VariableDescriptions{v},'(',''),')',''),'-','')
                  end
                  
                  
                  if b == 1
                      agg_table = t ;
                  else
                      agg_table = [agg_table ; t] ;
                  end
            end
            writetable(agg_table, fullfile( obj.meta.save_location, 'results', 'group_descriptives.xlsx')  );
            
            
            % Now build and write a wide format 
            
            
            
        end
        function simpleBoxPlots(obj)


%             pause() ;

            groupNum = obj.meta.PDVariableNames.groupNum;
            animal_id = obj.meta.PDVariableNames.animal_id;
            c_mean = obj.meta.PDVariableNames.mean;
            c_std = obj.meta.PDVariableNames.std;
            c_conmean = obj.meta.PDVariableNames.conmean;
            c_propmv = obj.meta.PDVariableNames.propmv;
            c_propNum = obj.meta.PDVariableNames.propNum;
            c_nrmMean = obj.meta.PDVariableNames.nrmMean;
            c_nrmConMean = obj.meta.PDVariableNames.nrmConMean;
            c_nrmPropMv = obj.meta.PDVariableNames.nrmPropMv;

            for b = 1 : numel(obj.processed.block)

                g = obj.returnGroups(obj.processed.block{b}.ind_desc(:,groupNum));
                figure;
                hold on;
                title(['mean ' obj.meta.dep_var],'Interpreter', 'none');
                boxplot(obj.processed.block{b}.ind_desc(:,c_mean),g)
                saveas(gcf,fullfile(obj.meta.save_location, 'results','oneway', obj.processed.block{b}.name, ['boxplot_mean.fig']));
                hold off;
                close(gcf);

            end



        end
        function R = zscore(obj,varargin)
            R = [] ;
            try
                %check inputs
                p = inputParser;
                p.addParameter('DV', ...
                    [],...
                    @(x) isnumeric(x));
                p.addParameter('GRP', ...
                    [],...
                    @(x) isnumeric(x));
                p.addParameter('PROP', ...
                    0,...
                    @(x) isnumeric(x));
                p.addParameter('SAV_LOC', ...
                    [],...
                    @(x) ischar(x));

                p.parse(varargin{:});
                inputs = p.Results;
            catch err
                obj.status = sprintf('FUNCTION "%s" CALL ERROR: %s', upper(obj.current_function),err.message);
                obj.help('onewayANOVA');
                R.error = 1;
                return;
            end


            if inputs.PROP == 1



            else
            end


        end
        function R = onewayANOVA(obj, varargin)
            R = [] ;
            try
                %check inputs
                p = inputParser;
                p.addParameter('DV', ...
                    [],...
                    @(x) isnumeric(x));
                p.addParameter('GRP', ...
                    [],...
                    @(x) isnumeric(x));
                p.addParameter('SAV_LOC', ...
                    [],...
                    @(x) ischar(x));

                p.parse(varargin{:});
                inputs = p.Results;
            catch err
                obj.status = sprintf('FUNCTION "%s" CALL ERROR: %s', upper(obj.current_function),err.message);
                obj.help('onewayANOVA');
                R.error = 1;
                return;
            end

            try
                [p,tbl,stats] = anova1(inputs.DV,inputs.GRP,'off');
                [c,m] = multcompare(stats,'Ctype','bonferroni','Display','off') ;
                % The first two columns of c show the groups that are compared. The fourth column shows the difference between the estimated group means. The third and fifth columns show the lower and upper limits for 95% confidence intervals for the true mean difference. The sixth column contains the p-value for a hypothesis test that the corresponding mean difference is equal to zero
                %             boxplot(obj.statistics.uni.BASELINE_mean,obj.statistics.uni.group_label2)

                ctab = array2table(c);
                ctab.Properties.VariableNames = {'Group1_num' 'Group2_num' 'LL' 'Difference' 'UL' 'p'} ;
                ctab.Group1 = stats.gnames(ctab.Group1_num);
                ctab.Group2 = stats.gnames(ctab.Group2_num);
                ctab = [ctab(:,7) ctab(:,8) ctab(:,1:6)];

                R.p = p ;
                R.tbl = tbl ;
                R.ctab = ctab ;

            catch err


                % doesn't work for one group, e.g.
                c = NaN;
                m = NaN;
                ctab = [];
            end


        end
        function sobj = saveobj(obj)
            % If the object does not have an account number,
            % Add account number to AccountNumber property

            % clear some unecessary data
            obj.data.merged_data = [];
            obj.data.viewpoint_data = [];
            obj.data.md_char = [] ;
            obj.data.md_num = [] ;
            obj.data.nbs = [] ;

            sobj = obj;
        end
        function listProcessed(obj)
            fprintf('NUM\tBLOCK\n') ;
            for o = 1 : numel(obj.processed.block)
                fprintf('%d\t%s\n',o,obj.processed.block{o}.name)
            end
        end
        function G= returnGroups(obj,grpNums)
            t = table(grpNums,'VariableNames', {'groupNum'});
                   C = join(t, obj.meta.group_cross);
            G = C.group_label2 ;
        end

    end

    methods (Access = protected)
        function splitData(obj)
            % Split into char and numeric
            % counts for initialization
            numCols = [] ;
            charCols = [] ;
            obj.meta.NVariableNames = struct('count', 0);
            obj.meta.CVariableNames = struct('count', 0);

            for w = 1 : width(obj.data.merged_data)
                thisCol = obj.data.merged_data.Properties.VariableNames{w} ;
                if isnumeric(obj.data.merged_data.(thisCol))
                    numCols = [numCols w];
                    %                     NVariableNames.count = NVariableNames.count + 1 ;
                    %                     NVariableNames.(obj.data.merged_data.Properties.VariableNames{w}) = NVariableNames.count ;
                    obj.registerNewDataColumn('NVariableNames', obj.data.merged_data.Properties.VariableNames{w});
                else
                    charCols = [charCols w] ;
                    %                    CVariableNames.count = NVariableNames.count + 1 ;
                    %                     CVariableNames.(obj.data.merged_data.Properties.VariableNames{w}) = CVariableNames.count ;
                    obj.registerNewDataColumn('CVariableNames', obj.data.merged_data.Properties.VariableNames{w});
                end

            end
            obj.data.md_num = obj.data.merged_data{:,numCols} ;
            obj.data.md_char = obj.data.merged_data{:,charCols} ;

            %             Build a group_num to group_name mapping
            obj.meta.group_cross = unique(obj.data.merged_data(:,{'groupNum' 'group_label2'}));



        end
        function labelEvents(obj, varargin)

            height = size(obj.data.md_num,1) ;
            xEnd = obj.meta.NVariableNames.xEnd ;
            start = obj.meta.NVariableNames.start ;


            NumEvents = sum(ismember(obj.data.timing_data.(obj.meta.event.event_index), obj.meta.event.event_label ));
            obj.meta.NumEvent = NumEvents ;
            Trunc_event_timing = obj.data.timing_data(ismember(obj.data.timing_data.(obj.meta.event.event_index), obj.meta.event.event_label ),:);

            event_iteration = repmat([0],height,1);
            event_preperiod = repmat([0],height,1);
            event_period = repmat([0],height,1);
            event_postperiod = repmat([0],height,1);


            for e = 1 : NumEvents
                preper_ons = Trunc_event_timing.Onset(e) - obj.meta.event.pre_event_period ;
                preper_off = Trunc_event_timing.Onset(e) ;

                per_ons = Trunc_event_timing.Onset(e) ;
                %                 per_off = Trunc_event_timing.Offset(e) + obj.meta.event.event_period ;

                per_off = Trunc_event_timing.Onset(e) + obj.meta.event.event_period ;

                postper_ons = per_off ;
                postper_off = per_off + obj.meta.event.post_event_period ;


                event_iteration(obj.data.md_num(:,xEnd) <=   postper_off &  obj.data.md_num(:,start) >=   preper_ons  ) = e ;
                event_preperiod(obj.data.md_num(:,xEnd) <=   preper_off &  obj.data.md_num(:,start) >=   preper_ons  )  = 1;
                event_period(obj.data.md_num(:,xEnd) <=   per_off &  obj.data.md_num(:,start) >=   per_ons  )  = 1;
                event_postperiod(obj.data.md_num(:,xEnd) <=   postper_off &  obj.data.md_num(:,start) >=   postper_ons  )   = 1;
            end

            obj.registerNewDataColumn('NVariableNames', 'event_iteration');
            obj.data.md_num=horzcat(obj.data.md_num,event_iteration);

            obj.registerNewDataColumn('NVariableNames', 'event_preperiod');
            obj.data.md_num=horzcat(obj.data.md_num,event_preperiod);

            obj.registerNewDataColumn('NVariableNames', 'event_period');
            obj.data.md_num=horzcat(obj.data.md_num,event_period);

            obj.registerNewDataColumn('NVariableNames', 'event_postperiod');
            obj.data.md_num=horzcat(obj.data.md_num,event_postperiod);

        end
        function restructureData(obj)
            % FUNCTION RESTRUCTUREDATA
            %   This restructures data into blocks defined by timing
            % ARGUMENTS: NONE
            %
            % NOTES:
            %



            %             reshape b blocks
            common_baseline = obj.meta.NVariableNames.common_baseline ;
            event_period = obj.meta.NVariableNames.event_period ;
            event_postperiod = obj.meta.NVariableNames.event_postperiod ;
            event_preperiod = obj.meta.NVariableNames.event_preperiod ;
            event_iteration = obj.meta.NVariableNames.event_iteration ;
            animal_id = obj.meta.NVariableNames.animal_id ;

            for event = 1 : obj.meta.NumEvent
                EVENT_DATA{event} = obj.data.md_num(obj.data.md_num(:,event_iteration)==event,:);
            end

            x = obj.data.md_num(:,common_baseline)==1 ;
            y = obj.data.md_num(:,event_period)==1 ;
            z = obj.data.md_num(:,event_postperiod)==1 ;

            block_data = [];
            block_data{1} = struct('name','ALL', 'data', obj.data.md_num(x | y | z,:)) ;
            block_data{2} = struct('name','BASE', 'data', obj.data.md_num(obj.data.md_num(:,common_baseline)==1,:)) ;
            block_data{3} = struct('name','STIM', 'data', obj.data.md_num(obj.data.md_num(:,event_period)==1,:)) ;
            block_data{4} = struct('name','POST', 'data', obj.data.md_num(obj.data.md_num(:,event_postperiod)==1,:)) ;
            for e = 1 : obj.meta.NumEvent
                block_data{end+1} = struct('name',['STIM' sprintf('%d',e)], 'data',...
                    obj.data.md_num((obj.data.md_num(:,event_period)==1) & (obj.data.md_num(:,event_iteration)==e),:  )) ;
                block_data{end+1} = struct('name',['POST' sprintf('%d',e)], 'data',...
                    obj.data.md_num((obj.data.md_num(:,event_postperiod)==1)& (obj.data.md_num(:,event_iteration)==e),:  )) ;
                block_data{end+1} = struct('name',['PRE' sprintf('%d',e)], 'data',...
                    obj.data.md_num((obj.data.md_num(:,event_preperiod)==1)& (obj.data.md_num(:,event_iteration)==e),:  )) ;
            end

            % stack by animals

            all_animals = unique(obj.data.md_num(:,animal_id)) ;
            nAnimals = numel(all_animals);



            nBlocks = numel(block_data) ;
            for b = 1 : nBlocks


                % determine largest stack size (hopefully all equal)
                maxSize = 0 ;
                for a = 1 : nAnimals
                    thisAnimal = all_animals(a) ;
                    thisSize = size(find(block_data{b}.data(:,animal_id)==thisAnimal),1) ;
                    if thisSize > maxSize
                        maxSize = thisSize ;
                    end
                end
                %                 pad if necessary
                thisWidth = size(block_data{b}.data(:,:),2);
                for a = 1 : nAnimals
                    thisAnimal = all_animals(a) ;
                    thisSize = size(find(block_data{b}.data(:,animal_id)==thisAnimal),1) ;
                    if thisSize < maxSize
                        discrepancy = maxSize -thisSize ;
                        pad = nan(discrepancy,thisWidth) ;
                        pad(:,animal_id) =  thisAnimal ;
                        block_data{b}.data=vertcat(block_data{b}.data,pad) ;

                    end

                end

                %               prealloc stacked data
                block_data{b}.stacked = nan(nAnimals,maxSize,thisWidth);

                %               Assign to stacked
                for a = 1 : nAnimals
                    thisAnimal = all_animals(a) ;
                    thisData =block_data{b}.data(block_data{b}.data(:,animal_id)==thisAnimal,:);
                    block_data{b}.stacked(a,:,:) = thisData;
                end
            end

            obj.data.nbs = block_data ;



        end
        function S = computeBasicStats(obj, varargin)
            % FUNCTION
            %
            % USAGE:
            % obj.computeBasicStats()
            %
            % ------------------------------------------------------------
            % REQUIRED ARGUMENTS:
            % 'ROWS',rows
            % 'DEP_VAR',dep_var
            %
            % OPTIONAL ARGUMENTS:
            %
            % OUTPUTS:
            %
            % NOTES:
            %


            try
                %check inputs
                p = inputParser;
                p.addParameter('ROWS', ...
                    [],...
                    @(x) isnumeric(x));
                p.addParameter('DEP_VAR', ...
                    [],...
                    @(x) ischar(x));
                p.addParameter('BY_EVENT', ...
                    0,...
                    @(x) isnumeric(x));

                %                 p.addParameter('DATA', ...
                %                     [],...
                %                     @(x) isstruct(x));
                p.parse(varargin{:});
                inputs = p.Results;
            catch err
                obj.status = sprintf('FUNCTION "%s" CALL ERROR: %s', upper(obj.current_function),err.message);
                obj.help('computeBasicStats');
                S.error = 1;
                return;
            end


            %subset table

            sTable = obj.data.merged_data(ismember(obj.data.merged_data.rowNum,inputs.ROWS),:) ;

            S = [];


            S.mean = mean(sTable.(inputs.DEP_VAR));
            S.std = std(sTable.(inputs.DEP_VAR));

            x = sTable.(inputs.DEP_VAR) ;
            x(x==0) = NaN ;
            S.conmean = nanmean(x) ;
            S.constd = nanstd(x) ;
            S.propMv = sum(sTable.(inputs.DEP_VAR)>0)/numel(sTable.(inputs.DEP_VAR)) ;
            S.propMvNum = sum(sTable.(inputs.DEP_VAR)>0);
            S.propMvDenom = numel(sTable.(inputs.DEP_VAR));

            thisAnimal = sTable.animal(1) ;
            priorRowNum = sTable.rowNum(1) - 1 ;
            if priorRowNum > 0
                priorRow = obj.data.merged_data(obj.data.merged_data.rowNum == priorRowNum,:) ;
                priorAnimal = priorRow.animal(1) ;
                if priorAnimal{1} == thisAnimal{1}
                    if priorRow.(inputs.DEP_VAR) > 0
                        S.premove = 1 ;
                    else
                        S.premove = 0 ;
                    end

                end
            else
                S.premove = NaN ;
            end

            if inputs.BY_EVENT
                %                 %if split by events, then
                %                 %determine the longest "event" array, probably all the
                %                 %same, but jsut double checking
                %                 maxL = 0;
                %                 for e = 1 : obj.meta.NumEvent
                %                     thisHeight = height(sTable(sTable.event_iteration==e,(inputs.DEP_VAR))) ;
                %                     if thisHeight > maxL
                %                         maxL = thisHeight ;
                %                     end
                %                 end
                %                 % create event arrays of length maxL
                %
                %                 event_arrays = zeros(maxL, obj.meta.NumEvent ) ;
                %
                %                 %populate
                %                 for e = 1 : obj.meta.NumEvent
                %                     thisArray = sTable{sTable.event_iteration==e,(inputs.DEP_VAR)} ;
                %                     for h = 1 : numel(thisArray)
                %                         event_arrays(h,e) = thisArray(h) ;
                %                     end
                %                 end
                %
                %                 S.ERA = mean(event_arrays,2);
                %

            else
                S.ERA = sTable.(inputs.DEP_VAR) ;
                S.bERA = sTable.(inputs.DEP_VAR)>0 ;
                S.cERA = S.ERA;
                S.cERA(S.bERA==0) = NaN ;
            end

            %             S.ERA_SM = movmean(S.ERA,[0 2]) ;


        end
        function C = computeTSparams(obj)


            obj.statistics.ts.BASELINE_mean = zeros(height(obj.statistics.ts),1) ;
            for h = 1 : height(obj.statistics.ts)
                obj.statistics.ts.BASELINE_mean(h) = mean(obj.statistics.ts{h,'BASELINE_ERA'}{1});

            end

            %             stitch together timing
            timeDim = 0 ;
            for h = 1 : height(obj.statistics.ts)
                td = numel(obj.statistics.ts.STIMPOST_S1_ERA{1}) + numel(obj.statistics.ts.STIMPOST_S2_ERA{1}) + numel(obj.statistics.ts.STIMPOST_S3_ERA{1}) + numel(obj.statistics.ts.STIMPOST_S4_ERA{1}) + numel(obj.statistics.ts.STIMPOST_S5_ERA{1}) ;
                if td > timeDim
                    timeDim = td ;
                end
            end
            %             for each group
            groups = unique(obj.statistics.ts.group_label2) ;

            numGroup = numel(groups) ;

            meanTimePlots = zeros(numGroup,timeDim) ;
            seTimePlots = zeros(numGroup,timeDim) ;

            eventBeats = zeros(1,5) ;

            % baselineTime = zeros(numGroup,1) ;

            for g = 1 : numGroup
                C(g).group = groups{g};
                temp = obj.statistics.ts(ismember(obj.statistics.ts.group_label2,groups{g}) ,:) ;
                bTime = zeros(height(temp),1) ;
                gTime = zeros(height(temp),timeDim) ;
                for a = 1 : height(temp)

                    eventBeats(1) =  0 ;
                    eventBeats(2) = numel(temp.STIMPOST_S1_ERA{a}) ;
                    eventBeats(3) = eventBeats(2) + numel(temp.STIMPOST_S2_ERA{a}) ;
                    eventBeats(4) = eventBeats(3) + numel(temp.STIMPOST_S3_ERA{a}) ;
                    eventBeats(5) = eventBeats(4) + numel(temp.STIMPOST_S4_ERA{a}) ;
                    eventBeats(6) = eventBeats(5) + numel(temp.STIMPOST_S5_ERA{a}) ;
                    gTime(a,:) = horzcat(temp.STIMPOST_S1_ERA{a},temp.STIMPOST_S2_ERA{a},temp.STIMPOST_S3_ERA{a},temp.STIMPOST_S4_ERA{a},temp.STIMPOST_S5_ERA{a});
                    bTime(a,:) = mean(temp.BASELINE_ERA{a}) ;
                end
                %                 meanTimePlots(g,:) = mean(gTime) - mean(bTime);
                meanTimePlots(g,:) = mean(gTime) ;

                %                 seTimePlots(g,:) = std(gTime)/(sqrt(height(temp) ));

                %                 L = Max
                %                 x0 = mid
                %                 k = steepness
                %                 c = y offset

                func = fittype('a*exp(-b*x)+c') ;
                %                 func = fittype('L*(1-1/(1+exp(-k*(x-x0))))+c');


                for seg = 1 : 5
                    %                     C.(groups{g}).(sprintf('seg%d',seg)) = 0 ;

                    %                   startPoints = [1,1,0];

                    segment = meanTimePlots(g,eventBeats(seg)+1: eventBeats(seg+1) )' ;
                    x = [0:numel(segment)-1]' ;
                    obj.meta.ts_seg = x ;
                    %param order = a,b,c
                    startPoints = [segment(1),...
                        1,...
                        segment(end),...
                        ];
                    %                    fprintf('%d\t%d\n',g,seg);
                    if g ==2 & seg == 4
                        %pause();
                    end
                    try
                        [fitobject,gof,output] = fit(x,segment,func,'Start', startPoints);
                    catch
                        fitobject.a = NaN;
                        fitobject.b = NaN;
                        fitobject.c = NaN;
                    end
                    %                     tempTime = [];
                    %                    for i = 1 : height(temp)
                    %                        tempTime = [tempTime;  obj.timeStitch(temp,i)] ;
                    %                    end

                    %                     upper = [200,10,20 ] ;
                    %                     lower = [] ;


                    %                     func = fittype('a*exp(-b*x)+c') ;
                    %                     [fitobject,gof,output] = fit(x,segment,func,'Start', startPoints,'Upper',upper ) ;
                    C(g).(sprintf('event_%d',seg)).a = fitobject.a ;
                    C(g).(sprintf('event_%d',seg)).b = fitobject.b ;
                    C(g).(sprintf('event_%d',seg)).c = fitobject.c ;
                    %                     figure ;
                    %                      hold on;
                    %                     p1 =  plot(x,segment) ;
                    %                     p2 = plot(fitobject) ;
                    %                     hold off ;
                    %
                    %                     xlim([0 numel(segment)-1] ) ;
                    % %                     p1.Color = rand(1,3) ;
                    %                     p2.Color = p2.rand(1,3) ;
                end
                %                 legend('1','2','3','4','5');
                %                 saveas(gcf,fullfile(obj.meta.save_location,'EventTimeSeries.fig'));
                %                 hold off ;
                %                 close(gcf);

            end

            obj.results.ts = C ;

            %%
            %             figure ;
            %             plot(segment) ;
            %             hold on ;





        end
        function T = printStatTable(obj)
            writetable(obj.statistics.uni, fullfile( obj.meta.save_location,'individualind_descs.xlsx')  );
            vn = obj.statistics.uni.Properties.VariableNames ;
            summ = grpstats(obj.statistics.uni,'group_label2',{'mean','sem'},'DataVars',   {    'BASELINE_mean'      ,...
                'BASELINE_std'       ,...
                'BASELINE_conmean'   ,...
                'BASELINE_propMv'    ,...
                'STIM_mean'          ,...
                'STIM_std'           ,...
                'STIM_conmean'       ,...
                'STIM_propMv'        ,...
                'POST_mean'          ,...
                'POST_std'           ,...
                'POST_conmean'       ,...
                'POST_propMv'        ,...
                'S1_mean'            ,...
                'S1_conmean'         ,...
                'S1_std'             ,...
                'S1_propMv'          ,...
                'S1_premove'         ,...
                'P1_mean'            ,...
                'P1_std'             ,...
                'P1_propMv'          ,...
                'P1_premove'         ,...
                'S2_mean'            ,...
                'S2_conmean'         ,...
                'S2_std'             ,...
                'S2_propMv'          ,...
                'S2_premove'         ,...
                'P2_mean'            ,...
                'P2_std'             ,...
                'P2_propMv'          ,...
                'P2_premove'         ,...
                'S3_mean'            ,...
                'S3_conmean'         ,...
                'S3_std'             ,...
                'S3_propMv'          ,...
                'S3_premove'         ,...
                'P3_mean'            ,...
                'P3_std'             ,...
                'P3_propMv'          ,...
                'P3_premove'         ,...
                'S4_mean'            ,...
                'S4_conmean'            ,...
                'S4_std'             ,...
                'S4_propMv'          ,...
                'S4_premove'         ,...
                'P4_mean'            ,...
                'P4_std'             ,...
                'P4_propMv'          ,...
                'P4_premove'         ,...
                'S5_mean'            ,...
                'S5_conmean'         ,...
                'S5_std'             ,...
                'S5_propMv'          ,...
                'S5_premove'         ,...
                'P5_mean'            ,...
                'P5_std'             ,...
                'P5_propMv'          ,...
                'P5_premove'         ,...
                'STIM_mean_nrm'      ,...
                'POST_mean_nrm'      ,...
                'POST_std_nrm'       ,...
                'S1_mean_nrm'        ,...
                'S2_mean_nrm'        ,...
                'S3_mean_nrm'        ,...
                'S4_mean_nrm'        ,...
                'S5_mean_nrm'        ,...
                'P1_mean_nrm'        ,...
                'P2_mean_nrm'        ,...
                'P3_mean_nrm'        ,...
                'P4_mean_nrm'        ,...
                'P5_mean_nrm'        ,...
                'BASELINEvsSTIM_mean',...
                'STIMvPOST_mean'     ,...
                'S1vS5_mean'         } ) ;
            writetable(summ, fullfile( obj.meta.save_location,'groupind_descs.xlsx')  );
            obj.basicTimeSeriesPlots();
            %             obj.basicProporionPlots() ;
            ptab = table;

            fn = fieldnames(obj.results.oneway);
            %             for f = 1 : numel(fn)
            %                 obj.results.oneway.(fn{f}).type
            %             end
            for f = 1 : numel(fn)
                if strcmp(obj.results.oneway.(fn{f}).type,'onewayANOVA')
                    try
                        rtab = cell2table(obj.results.oneway.(fn{f}).tbl(2:end,2:end));
                    catch
                    end

                    rtab.Properties.VariableNames = genvarname(obj.results.oneway.(fn{f}).tbl(1,2:end)) ;
                    rtab.Properties.VariableNames{end} = 'P' ;
                    rtab.Properties.RowNames = genvarname(obj.results.oneway.(fn{f}).tbl(2:end,1)');
                    writetable(rtab, fullfile( obj.meta.save_location,'oneway',(fn{f}),'onewayTable.xlsx' ) ,'WriteRowNames' , 1);
                    newRow = cell2table({'oneway', {fn{f}}, rtab.P('Groups')},'VariableNames',{'type','var','p_group'});
                    ptab = [ptab;newRow];
                    if ~isempty(obj.results.oneway.(fn{f}).mult_c)
                        writetable(obj.results.oneway.(fn{f}).mult_c, fullfile( obj.meta.save_location,'oneway',(fn{f}),'multipleComparisons.xlsx' ) ,'WriteRowNames' , 1);
                    end
                end
            end
            writetable(ptab, fullfile( obj.meta.save_location,'pSummaryOneway.xlsx' ) );


            ptab = table;
            fn = fieldnames(obj.results.twoway);
            for f = 1 : numel(fn)
                writetable(obj.results.twoway.(fn{f}).btw_tbl,fullfile( obj.meta.save_location,'twoway',(fn{f}),'betweenTable.xlsx') ,'WriteRowNames' , 1) ;
                try
                    writetable(obj.results.twoway.(fn{f}).within_tbl,fullfile( obj.meta.save_location,'twoway',(fn{f}),'withinTable.xlsx') ,'WriteRowNames' , 1 );

                    wpval = obj.results.twoway.(fn{f}).within_tbl{'Group:Within','pValue'} ;
                catch
                    % account for missing within_table due to single group
                    wpval = nan;
                end
                bpval = obj.results.twoway.(fn{f}).btw_tbl.pValue(2) ;

                newRow = cell2table({'twoway', {fn{f}}, bpval,wpval},'VariableNames',{'type','var','p_btw','p_interaction'});
                ptab = [ptab;newRow];

            end
            writetable(ptab, fullfile( obj.meta.save_location,'pSummaryTwoway.xlsx' ) );
            %             1.	STIM_mean: 1 way
            %             2.	POST_mean: 1 way
            %             3.	S1_mean: 1 way
            %             4.	S2_mean: 1 way
            %             5.	S3_mean: 1 way
            %             6.	S4_mean: 1 way
            %             7.	S5_mean: 1 way
            %             8.	P1_std: 1 way
            %             9.	P2_std: 1 way
            %             10.	P3_std: 1 way
            %             11.	P4_std: 1 way
            %             12.	P5_std: 1 way
            %             13.	STIM-POST : 2 way, derived
            %             14.	S1-S5_mean : 2 way, derived
            %             15.	P1-P5_mean : 2 way, derived


            % make figures for ts estimates

            groups = {obj.results.ts.group};

            %             for e = 1 : obj.meta.NumEvent
            %                 figure;
            %                 hold on;
            %                 for f = 1 : numel(obj.results.ts)
            %                     a = obj.results.ts(f).(sprintf('event_%d',e)).a ;
            %                     b = obj.results.ts(f).(sprintf('event_%d',e)).b ;
            %                     c = obj.results.ts(f).(sprintf('event_%d',e)).c ;
            %                     x = obj.meta.ts_seg;
            %                     y = a*exp(-b*x)+c ;
            %                     plot(x,y)
            %
            %                 end
            %                 legend(groups)
            %                 saveas(gcf,fullfile(obj.meta.save_location,sprintf('fittedEventPlot_%d.fig',e)));
            %                 hold off;
            %                 close(gcf);
            %             end

        end


        function R = onewayANOVA2(obj, varargin)
            % FUNCTION
            %  Test the significance of group on single univariate stat
            %  (e.g do groups differ on mean motion during post-event response)
            % USAGE:
            % obj.computeBasicStats()
            %
            % ------------------------------------------------------------
            % REQUIRED ARGUMENTS:
            % 'ROWS',rows
            % 'DEP_VAR',dep_var
            %
            % OPTIONAL ARGUMENTS:
            %
            % OUTPUTS:
            %
            % NOTES:
            %

            R = [] ;

            try
                %check inputs
                p = inputParser;
                p.addParameter('VAR', ...
                    [],...
                    @(x) ischar(x));
                p.addParameter('GRP_COL', ...
                    [],...
                    @(x) ischar(x));
                p.addParameter('TITLE', ...
                    [],...
                    @(x) ischar(x));
                p.addParameter('NormTo', ...
                    [],...
                    @(x) ischar(x));
                p.parse(varargin{:});
                inputs = p.Results;
            catch err
                obj.status = sprintf('FUNCTION "%s" CALL ERROR: %s', upper(obj.current_function),err.message);
                obj.help('onewayANOVA');
                R.error = 1;
                return;
            end

            mkdir(fullfile(obj.meta.save_location,'oneway')) ;
            mkdir(fullfile(obj.meta.save_location,'oneway',inputs.TITLE)) ;
            if ~isempty(inputs.NormTo)
                try
                    [p,tbl,stats] = anova1(obj.statistics.uni.(inputs.VAR) -obj.statistics.uni.(inputs.NormTo),obj.statistics.uni.(inputs.GRP_COL),'off');
                catch err
                    obj.status = sprintf('FUNCTION "%s" NORMTO ERROR: %s', upper(obj.current_function),err.message);
                    return ;
                end
            else
                [p,tbl,stats] = anova1(obj.statistics.uni.(inputs.VAR),obj.statistics.uni.(inputs.GRP_COL),'off');
            end

            try
                [c,m] = multcompare(stats,'Ctype','bonferroni','Display','off') ;
                % The first two columns of c show the groups that are compared. The fourth column shows the difference between the estimated group means. The third and fifth columns show the lower and upper limits for 95% confidence intervals for the true mean difference. The sixth column contains the p-value for a hypothesis test that the corresponding mean difference is equal to zero
                %             boxplot(obj.statistics.uni.BASELINE_mean,obj.statistics.uni.group_label2)

                ctab = array2table(c);
                ctab.Properties.VariableNames = {'Group1_num' 'Group2_num' 'LL' 'Difference' 'UL' 'p'} ;
                ctab.Group1 = stats.gnames(ctab.Group1_num);
                ctab.Group2 = stats.gnames(ctab.Group2_num);
                ctab = [ctab(:,7) ctab(:,8) ctab(:,1:6)];
            catch


                % doesn't work for one group, e.g.
                c = NaN;
                m = NaN;
                ctab = [];
            end



            figure;
            hold on;
            title(inputs.VAR,'Interpreter', 'none');
            boxplot(obj.statistics.uni.(inputs.VAR),obj.statistics.uni.group_label2);
            saveas(gcf,fullfile(obj.meta.save_location, 'oneway',inputs.TITLE, [inputs.VAR '.fig']));
            hold off;
            close(gcf);






            R.type = 'onewayANOVA' ;
            R.var = inputs.VAR;
            R.p = p;
            R.tbl = tbl;
            R.stats = stats ;
            R.mult_c = ctab ;
            R.boxplot  = fullfile(obj.meta.save_location, 'oneway', inputs.TITLE, [inputs.VAR '.fig']) ;


            obj.results.oneway.(genvarname(inputs.TITLE)) = R ;

        end
        function R = onewayBAR(obj, varargin)
            % FUNCTION
            %  Test the significance of group on single univariate stat
            %  (e.g do groups differ on mean motion during post-event response)
            % USAGE:
            % obj.computeBasicStats()
            %
            % ------------------------------------------------------------
            % REQUIRED ARGUMENTS:
            % 'ROWS',rows
            % 'DEP_VAR',dep_var
            %
            % OPTIONAL ARGUMENTS:
            %
            % OUTPUTS:
            %
            % NOTES:
            %

            R = [] ;

            try
                %check inputs
                p = inputParser;
                p.addParameter('VAR', ...
                    [],...
                    @(x) ischar(x));
                p.addParameter('GRP_COL', ...
                    [],...
                    @(x) ischar(x));
                p.addParameter('TITLE', ...
                    [],...
                    @(x) ischar(x));
                p.addParameter('NormTo', ...
                    [],...
                    @(x) ischar(x));
                p.parse(varargin{:});
                inputs = p.Results;
            catch err
                obj.status = sprintf('FUNCTION "%s" CALL ERROR: %s', upper(obj.current_function),err.message);
                obj.help('onewayBAR');
                R.error = 1;
                return;
            end

            mkdir(fullfile(obj.meta.save_location,'oneway')) ;
            mkdir(fullfile(obj.meta.save_location,'oneway',inputs.TITLE)) ;


            summ = grpstats(obj.statistics.uni,inputs.GRP_COL,{'mean','sem'},'DataVars',inputs.VAR);
            figure;
            hold on;
            errorbar(1:numel(summ.(['mean_' inputs.VAR])),...
                summ.(['mean_' inputs.VAR]),...
                summ.(['sem_' inputs.VAR]),...
                'd') ;
            xlim([0,numel(summ.(['mean_' inputs.VAR]))+1]);
            xticks([0:numel(summ.(['mean_' inputs.VAR]))+1]);

            xticklabels( horzcat({''},summ.group_label2',{''}));
            title(inputs.VAR,'Interpreter', 'none');
            saveas(gcf,fullfile(obj.meta.save_location, 'oneway',inputs.TITLE, [inputs.VAR '.fig']));
            hold off;
            close(gcf);




            R.type = 'onewayBAR' ;
            R.var = inputs.VAR;
            R.barplot  = fullfile(obj.meta.save_location, 'oneway', inputs.TITLE, [inputs.VAR '.fig']) ;


            obj.results.oneway.(genvarname(inputs.TITLE)) = R ;

        end
        function R = onewayCHI(obj, varargin)
            % FUNCTION
            %  Test the significance of group on single univariate stat
            %  (e.g do groups differ on mean motion during post-event response)
            % USAGE:
            % obj.computeBasicStats()
            %
            % ------------------------------------------------------------
            % REQUIRED ARGUMENTS:
            % 'ROWS',rows
            % 'DEP_VAR',dep_var
            %
            % OPTIONAL ARGUMENTS:
            %
            % OUTPUTS:
            %
            % NOTES:
            %

            R = [] ;

            try
                %check inputs
                p = inputParser;
                p.addParameter('VAR', ...
                    [],...
                    @(x) ischar(x));
                p.addParameter('GRP_COL', ...
                    [],...
                    @(x) ischar(x));
                p.addParameter('TITLE', ...
                    [],...
                    @(x) ischar(x));
                p.addParameter('NormTo', ...
                    [],...
                    @(x) ischar(x));
                p.parse(varargin{:});
                inputs = p.Results;
            catch err
                obj.status = sprintf('FUNCTION "%s" CALL ERROR: %s', upper(obj.current_function),err.message);
                obj.help('onewaychi');
                R.error = 1;
                return;
            end


            mkdir(fullfile(obj.meta.save_location,'oneway')) ;
            mkdir(fullfile(obj.meta.save_location,'oneway',inputs.TITLE)) ;


            summ = grpstats(obj.statistics.uni,inputs.GRP_COL,{'mean','sem'},'DataVars',inputs.VAR);
            figure;
            hold on;
            errorbar(1:numel(summ.(['mean_' inputs.VAR])),...
                summ.(['mean_' inputs.VAR]),...
                summ.(['sem_' inputs.VAR]),...
                'd') ;
            xlim([0,numel(summ.(['mean_' inputs.VAR]))+1]);
            xticks([0:numel(summ.(['mean_' inputs.VAR]))+1]);

            xticklabels( horzcat({''},summ.group_label2',{''}));
            title(inputs.VAR,'Interpreter', 'none');
            saveas(gcf,fullfile(obj.meta.save_location, 'oneway',inputs.TITLE, [inputs.VAR '.fig']));
            hold off;
            close(gcf);




            R.type = 'onewayCHI' ;
            R.var = inputs.VAR;
            R.barplot  = fullfile(obj.meta.save_location, 'oneway', inputs.TITLE, [inputs.VAR '.fig']) ;


            obj.results.oneway.(genvarname(inputs.TITLE)) = R ;




        end
        function R =  twowayANOVA(obj,varargin )
            % FUNCTION
            %   Test the significance of group (between subjects), event components (within subjects), and their
            %   interaction  on single univariate stat (e.g does mean
            %   motion differe between groups, between event and post-event, or the interaction )
            % NOTE: the within subjects variables must be specified as
            % separate columns in the statistics.uni table


            R = [] ;

            try
                %check inputs
                p = inputParser;
                p.addParameter('VAR1', ...
                    [],...
                    @(x) ischar(x));

                p.addParameter('VAR2', ...
                    [],...
                    @(x) ischar(x));
                p.addParameter('GRP_COL', ...
                    [],...
                    @(x) ischar(x));
                p.addParameter('TITLE', ...
                    [],...
                    @(x) ischar(x));
                p.addParameter('NormTo', ...
                    [],...
                    @(x) ischar(x));
                p.parse(varargin{:});
                inputs = p.Results;
            catch err
                obj.status = sprintf('FUNCTION "%s" CALL ERROR: %s', upper(obj.current_function),err.message);
                obj.help('twowayANOVA');
                R.error = 1;
                return;
            end


            mkdir(fullfile(obj.meta.save_location,'twoway') );
            mkdir(fullfile(obj.meta.save_location,'twoway',inputs.TITLE)) ;
            % Construct rm table spec %
            rmtab = obj.statistics.uni(:,{'group_label2',inputs.VAR1,inputs.VAR2}) ;

            if ~isempty(inputs.NormTo)
                try
                    rmtab.(inputs.VAR1) = rmtab.(inputs.VAR1) - obj.statistics.uni.(inputs.NormTo);
                    rmtab.(inputs.VAR2) = rmtab.(inputs.VAR2) - obj.statistics.uni.(inputs.NormTo);
                catch err
                    obj.status = sprintf('FUNCTION "%s" NORMTO ERROR: %s', upper(obj.current_function),err.message);
                    return ;
                end
            end
            rmtab.Properties.VariableNames = {'Group', 'y1','y2'} ;
            % construct within design spec %

            rm = fitrm(rmtab,'y1-y2~Group','WithinDesign',[1 2]');


            %within subj tests
            try
                rt = rm.ranova();

                rt.Properties.RowNames = {'(Intercept):Within' ; 'Group:Within' ; 'Error(Within)' } ;
            catch
                % one group?
                rt = [];
            end

            %between subj tests

            bt = rm.anova() ;


            %box plot
            elab = [repmat({inputs.VAR1},height(rmtab),1); repmat({inputs.VAR2},height(rmtab),1)];
            glab = repmat(rmtab.Group,2,1) ;


            figure;
            hold on;
            title([inputs.VAR1 ' vs ' inputs.VAR2] ,'Interpreter', 'none');
            boxplot([rmtab.y1; rmtab.y2],{glab, elab}, 'factorgap',10, 'LabelOrientation','inline')
            saveas(gcf,fullfile(obj.meta.save_location, 'twoway',inputs.TITLE, [inputs.VAR1 'vs' inputs.VAR2 '.fig']));
            hold off;
            close(gcf);


            %             % Suppose we want to compare levels of Attention for each combination
            %             % of levels of TestCond and TMS.
            %
            %             % 1. Convert factors to categorical.
            %             within2 = within;
            %             within2.Attention = categorical(within2.Attention);
            %             within2.TestCond = categorical(within2.TestCond);
            %             within2.TMS = categorical(within2.TMS);
            %
            %             % 2. Create an interaction factor capturing each combination of levels
            %             % of TestCond and TMS.
            %             within2.TestCond_TMS = within2.TestCond .* within2.TMS;
            %
            %             % 3. Call fitrm with the modified within design.
            %             rm2 = fitrm(t,'Y1-Y8~1','WithinDesign',within2);
            %             ranovatbl2 = ranova(rm2, 'WithinModel','TestCond*Attention*TMS')
            %
            %             % 4. Use interaction factor TestCond_TMS as the 'By' variable in multcompare.
            %             multcompare(rm2,'Attention','By','TestCond_TMS')


            % mult comparison
            mtab = multcompare(rm,'Time','By','Group') ;
            mtab.Properties.VariableNames = {'Group'    inputs.VAR1    inputs.VAR2    'Difference'    'StdErr'      'pValue'       'Lower'      'Upper' };


            R.type = 'twowayANOVA' ;
            R.var1 = inputs.VAR1;
            R.var2 = inputs.VAR2;
            R.within_tbl = rt;
            R.btw_tbl = bt ;
            %             R.stats = stats ;

            R.mult_c =  mtab ;
            R.boxplot  = fullfile(obj.meta.save_location, 'twoway', inputs.TITLE, [inputs.VAR1 'vs' inputs.VAR2 '.fig']) ;

            t = R.mult_c(2:2:end,:);
            figure ;
            hold on;
            x = 1:numel(t.Difference) ;
            bar(x,t.Difference);
            xticks(x) ;
            xticklabels(t.Group');
            errorbar(x,t.Difference,t.StdErr,'.');
            saveas(gcf,fullfile(obj.meta.save_location, 'twoway',inputs.TITLE, [inputs.VAR1 'vs' inputs.VAR2 '_diff.fig']));
            R.diffPlot = fullfile(obj.meta.save_location, 'twoway',inputs.TITLE, [inputs.VAR1 'vs' inputs.VAR2 '_diff.fig']) ;

            hold off;
            close(gcf);
            %
            obj.results.twoway.(genvarname(inputs.TITLE)) = R ;

        end
        function p = basicTimeSeriesPlots(obj,varargin)
            try
                %check inputs
                p = inputParser;

                p.addParameter('TITLE', ...
                    [],...
                    @(x) ischar(x));

                p.parse(varargin{:});
                inputs = p.Results;
            catch err
                obj.status = sprintf('FUNCTION "%s" CALL ERROR: %s', upper(obj.current_function),err.message);
                obj.help('basicTimeSeriesPlots');
                p.error = 1;
                return;
            end

            %define groups
            groups = unique(obj.statistics.ts.group_label2) ;
            numGroup = numel(groups) ;

            timeDim = 0 ;
            for h = 1 : height(obj.statistics.ts)
                td = numel(obj.statistics.ts.BASELINE_ERA{h})+ numel(obj.statistics.ts.STIMPOST_S1_ERA{h}) + numel(obj.statistics.ts.STIMPOST_S2_ERA{h}) + numel(obj.statistics.ts.STIMPOST_S3_ERA{h}) + numel(obj.statistics.ts.STIMPOST_S4_ERA{h}) + numel(obj.statistics.ts.STIMPOST_S5_ERA{h}) ;
                if td > timeDim
                    timeDim = td ;
                end
            end

            for g = 1 : numGroup
                temp = obj.statistics.ts(ismember(obj.statistics.ts.group_label2,groups{g}) ,:) ;
                %                 bTime = zeros(height(temp),1) ;
                %                 structures for aggregrating
                gTime = zeros(height(temp),timeDim) ;
                bgTime = zeros(height(temp),timeDim) ;
                congTime = zeros(height(temp),timeDim) ;
                for a = 1 : height(temp)
                    gTime(a,:) = horzcat(temp.BASELINE_ERA{a},temp.STIMPOST_S1_ERA{a},temp.STIMPOST_S2_ERA{a},temp.STIMPOST_S3_ERA{a},temp.STIMPOST_S4_ERA{a},temp.STIMPOST_S5_ERA{a});
                    bgTime(a,:) = horzcat(temp.BASELINE_bERA{a},temp.STIMPOST_S1_bERA{a},temp.STIMPOST_S2_bERA{a},temp.STIMPOST_S3_bERA{a},temp.STIMPOST_S4_bERA{a},temp.STIMPOST_S5_bERA{a});
                    congTime(a,:) = horzcat(temp.BASELINE_cERA{a},temp.STIMPOST_S1_cERA{a},temp.STIMPOST_S2_cERA{a},temp.STIMPOST_S3_cERA{a},temp.STIMPOST_S4_cERA{a},temp.STIMPOST_S5_cERA{a});
                end
                meanTimePlots(g,:) = mean(gTime,1) ;
                bmeanTimePlots(g,:) = mean(bgTime,1) ;
                conmeanTimePlots(g,:) = nanmean(congTime,1) ;

                seTimePlots(g,:) = std(gTime,0,1)/(sqrt(height(temp) ));
                bseTimePlots(g,:) = zeros(1,size(gTime,2)) ;
                conseTimePlots(g,:) = nanstd(congTime,0,1)/(sqrt(height(temp) ));
            end


            % FIGUREWITHBASELINE   NOT NORMD based on mean activity

            figure;
            hold on;
            title(inputs.TITLE,'Interpreter', 'none');
            cmap = lines(numGroup);
            for g = numGroup:-1:1
                %compute timeise SEM per group
                y = meanTimePlots(g,:);
                x = 1:numel(y);
                sem = seTimePlots(g,:) ;
                boundedline(x,y,sem,'alpha','cmap', cmap(g,:), 'transparency', 0.5);
            end
            h = get(gca,'Children');
            h1 = h(1:2:end);
            legend(h1,groups,'Interpreter','none');

            saveas(gcf,fullfile(obj.meta.save_location,'EventTimeSeries_withbaseline_raw.fig'));
            hold off;
            close(gcf);

            % FIGUREWITHBASELINE   NOT NORMD based on binary activity
            figure;
            hold on;
            title(inputs.TITLE,'Interpreter', 'none');
            cmap = lines(numGroup);
            for g = numGroup:-1:1
                %compute timeise SEM per group
                y = bmeanTimePlots(g,:);
                x = 1:numel(y);
                sem = bseTimePlots(g,:) ;
                %                 boundedline(x,y,sem,'alpha','cmap', cmap(g,:), 'transparency', 0.5);
                subplot(numGroup,1,g) ;
                plot(x,y);
                legend(groups{g},'location','northwest');
            end

            saveas(gcf,fullfile(obj.meta.save_location,'BinaryEventTimeSeries_withbaseline.fig'));
            hold off;
            close(gcf);

            %FIGUREWITHBASELINE   NOT NORMD based on contingent activity
            figure;
            hold on;
            title(inputs.TITLE,'Interpreter', 'none');
            cmap = lines(numGroup);
            for g = numGroup:-1:1
                %compute timeise SEM per group
                y = conmeanTimePlots(g,:);
                x = 1:numel(y);
                sem = conseTimePlots(g,:) ;
                %                 boundedline(x,y,sem,'alpha','cmap', cmap(g,:), 'transparency', 0.5);
                subplot(numGroup,1,g) ;
                plot(x,y);
                legend(groups{g},'location','northwest');
            end

            saveas(gcf,fullfile(obj.meta.save_location,'contingentEventTimeSeries_withbaseline.fig'));
            hold off;
            close(gcf);



            %   EVENT ONLY FIGURE, NORMED

            %             compute baseline mean values

            obj.statistics.ts.BASELINE_mean = zeros(height(obj.statistics.ts),1) ;
            for h = 1 : height(obj.statistics.ts)
                obj.statistics.ts.BASELINE_mean(h) = mean(obj.statistics.ts{h,'BASELINE_ERA'}{1});

            end

            %             stitch together timing
            timeDim = 0 ;
            for h = 1 : height(obj.statistics.ts)
                td = numel(obj.statistics.ts.STIMPOST_S1_ERA{h}) + numel(obj.statistics.ts.STIMPOST_S2_ERA{h}) + numel(obj.statistics.ts.STIMPOST_S3_ERA{h}) + numel(obj.statistics.ts.STIMPOST_S4_ERA{h}) + numel(obj.statistics.ts.STIMPOST_S5_ERA{h}) ;
                if td > timeDim
                    timeDim = td ;
                end
            end

            %             Reset structures
            meanTimePlots = zeros(numGroup,timeDim) ;
            bmeanTimePlots = zeros(numGroup,timeDim) ;

            seTimePlots = zeros(numGroup,timeDim) ;

            eventBeats = zeros(1,5) ;

            % baselineTime = zeros(numGroup,1) ;

            for g = 1 : numGroup
                temp = obj.statistics.ts(ismember(obj.statistics.ts.group_label2,groups{g}) ,:) ;
                bTime = zeros(height(temp),1) ;
                gTime = zeros(height(temp),timeDim) ;
                binTime = zeros(height(temp),timeDim) ;

                for a = 1 : height(temp)
                    eventBeats(1) =  0 ;
                    eventBeats(2) = numel(temp.STIMPOST_S1_ERA{a}) ;
                    eventBeats(3) = eventBeats(2) + numel(temp.STIMPOST_S2_ERA{a}) ;
                    eventBeats(4) = eventBeats(3) + numel(temp.STIMPOST_S3_ERA{a}) ;
                    eventBeats(5) = eventBeats(4) + numel(temp.STIMPOST_S4_ERA{a}) ;
                    gTime(a,:) = horzcat(temp.STIMPOST_S1_ERA{a},temp.STIMPOST_S2_ERA{a},temp.STIMPOST_S3_ERA{a},temp.STIMPOST_S4_ERA{a},temp.STIMPOST_S5_ERA{a});
                    binTime(a,:) = horzcat(temp.STIMPOST_S1_bERA{a},temp.STIMPOST_S2_bERA{a},temp.STIMPOST_S3_bERA{a},temp.STIMPOST_S4_bERA{a},temp.STIMPOST_S5_bERA{a});
                    bTime(a,:) = mean(temp.BASELINE_ERA{a}) ;
                end
                normed_meanTimePlots(g,:) = mean(gTime) - mean(bTime);
                unnormed_meanTimePlots(g,:) = mean(gTime) ;
                bmeanTimePlots(g,:) = mean(binTime) ;
                seTimePlots(g,:) = std(gTime)/(sqrt(height(temp) ));
            end


            %             plotting

            groups = unique(obj.statistics.ts.group_label2) ;
            numGroup = numel(groups) ;

            %             Event period only basline normed

            figure;
            hold on;
            title(inputs.TITLE,'Interpreter', 'none');
            cmap = lines(numGroup);
            for g = numGroup:-1:1
                %compute timeise SEM per group
                y = normed_meanTimePlots(g,:);
                x = 1:numel(y);
                sem = seTimePlots(g,:) ;
                boundedline(x,y,sem,'alpha','cmap', cmap(g,:), 'transparency', 0.5);
            end

            %Event Raster
            eFlag = 0;
            yl = ylim;
            raster_xs = eventBeats;
            raster_ys =  repmat(yl(1),size(raster_xs,1),1);
            plot(raster_xs, raster_ys,'ro');
            h = get(gca,'Children');

            legend([h(6:2:end)],[groups],'Interpreter','none');

            saveas(gcf,fullfile(obj.meta.save_location,'EventTimeSeries_normed.fig'));
            hold off;
            close(gcf);



            %             Event period only unnormed normed

            figure;
            hold on;
            title(inputs.TITLE,'Interpreter', 'none');
            cmap = lines(numGroup);
            for g = numGroup:-1:1
                %compute timeise SEM per group
                y = unnormed_meanTimePlots(g,:);
                x = 1:numel(y);
                sem = seTimePlots(g,:) ;
                boundedline(x,y,sem,'alpha','cmap', cmap(g,:), 'transparency', 0.5);
            end

            %Event Raster
            eFlag = 0;
            yl = ylim;
            raster_xs = eventBeats;
            raster_ys =  repmat(yl(1),size(raster_xs,1),1);
            plot(raster_xs, raster_ys,'ro');
            h = get(gca,'Children');

            legend([h(6:2:end)],[groups],'Interpreter','none');
            saveas(gcf,fullfile(obj.meta.save_location,'EventTimeSeries_raw.fig'));
            hold off;
            close(gcf);


            %           Binary activity event period only


            figure;
            hold on;
            title(inputs.TITLE,'Interpreter', 'none');
            cmap = lines(numGroup);
            for g = numGroup:-1:1
                %compute timeise SEM per group
                y = bmeanTimePlots(g,:);
                x = 1:numel(y);
                sem = seTimePlots(g,:) ;
                subplot(numGroup,1,g) ;
                hold on ;
                plot(x,y);
                %Event Raster
                eFlag = 0;
                yl = ylim;
                raster_xs = eventBeats;
                raster_ys =  repmat(yl(1),size(raster_xs,1),1);
                plot(raster_xs, raster_ys,'ro');
                legend(groups{g},'location','eastoutside');
                hold off ;
            end

            saveas(gcf,fullfile(obj.meta.save_location,'BinaryEventTimeSeries.fig'));
            hold off;
            close(gcf);

        end

        function computeDerivedVariables(obj)

            animal_id = obj.meta.PDVariableNames.animal_id;
            groupNum =obj.meta.PDVariableNames.groupNum;
            cmean = obj.meta.PDVariableNames.mean;
            cstd = obj.meta.PDVariableNames.std;
            cconmean = obj.meta.PDVariableNames.conmean;
            cpropmv = obj.meta.PDVariableNames.propmv;
            cpropNum = obj.meta.PDVariableNames.propNum;


            %             register new variables
            obj.registerNewDataColumn('PDVariableNames', 'nrmMean');
            obj.registerNewDataColumn('PDVariableNames', 'nrmConMean');
            obj.registerNewDataColumn('PDVariableNames', 'nrmPropMv');
            %             we compute baseline norming for all blocks, including
            %             baseline so there there is consistency across the results
            %             matrices

            %             Add new blocks defined by specific comparisons
            %             S1 v S5
            obj.processed.block{end+1} = struct('name','S1-S5',...
                'ind_desc', obj.processed.block{5}.ind_desc    - obj.processed.block{17}.ind_desc,...
                'ind_timeseries',obj.processed.block{5}.ind_timeseries    - mean(obj.processed.block{17}.ind_timeseries(:,3:end),2));
            % fill in ids and grp numbers
            obj.processed.block{end}.ind_desc(:,1:2) = obj.processed.block{5}.ind_desc(:,1:2);
            obj.processed.block{end}.ind_timeseries(:,1:2) = obj.processed.block{5}.ind_timeseries(:,1:2);


            obj.processed.block{end+1} = struct('name','STIM-POST',...
                'ind_desc', obj.processed.block{3}.ind_desc    - obj.processed.block{4}.ind_desc,...
                'ind_timeseries',obj.processed.block{3}.ind_timeseries    - mean(obj.processed.block{4}.ind_timeseries(:,3:end),2) );
            obj.processed.block{end}.ind_desc(:,1:2) = obj.processed.block{3}.ind_desc(:,1:2);
            obj.processed.block{end}.ind_timeseries(:,1:2) = obj.processed.block{2}.ind_timeseries(:,1:2);

            for b = 1 : numel(obj.processed.block)
                %                 Compute new variables

                nrmMean = obj.processed.block{b}.ind_desc(:,cmean) -  obj.processed.block{2}.ind_desc(:,cmean);
                nrmConMean = obj.processed.block{b}.ind_desc(:,cconmean) -  obj.processed.block{2}.ind_desc(:,cconmean);
                nrmPropMv = obj.processed.block{b}.ind_desc(:,cpropmv) -  obj.processed.block{2}.ind_desc(:,cpropmv);
                %                 nrmPropNum= obj.processed.block{b}.data(:,cmean) -  obj.processed.block{1}.data(:,cmean);
                obj.processed.block{b}.ind_desc = [obj.processed.block{b}.ind_desc nrmMean nrmConMean nrmPropMv];

            end

        end
        function loadStruct = loadFiles(obj,varargin)
            % FUNCTION
            %
            % USAGE:
            % obj.function()
            %
            % ------------------------------------------------------------
            % REQUIRED ARGUMENTS:
            % 'HANDLE_PARENT',parent
            % 'KEY_ID',key_id)
            %
            % OPTIONAL ARGUMENTS:
            %
            % OUTPUTS:
            %
            % NOTES:
            %

            loadStruct = struct('error', 0);

            validFileSet = {'GROUPING' 'VIEWPOINT' 'TIMING'};
            try
                %check inputs
                p = inputParser;
                p.addParameter('FILE_SET', ...
                    [],...
                    @(x) ismember(x,validFileSet));
                p.parse(varargin{:});
                inputs = p.Results;
            catch err
                obj.status = sprintf('FUNCTION "%s" CALL ERROR: %s', upper(obj.current_function),err.message);
                obj.help('loadFiles');
                return;
            end

            switch inputs.FILE_SET
                case 'GROUPING'
                    loadFunc = 'loadGroupingFile';
                    loadFile = obj.files.grouping_file;
                    loadTable = 'group_data';
                case 'VIEWPOINT'
                    loadFunc = 'loadViewpointFile';
                    loadFile = obj.files.viewpoint_file;
                    loadTable = 'viewpoint_data';
                case 'TIMING'
                    loadFunc = 'loadTimingFile';
                    loadFile = obj.files.timing_file;
                    loadTable = 'timing_data';
            end


            if ischar(loadFile)  % definitely only 1 file
                dStruct = obj.(loadFunc)('FILE',loadFile);
                if dStruct.error
                    loadStruct.error = 1;
                    return;
                end
                obj.data.(loadTable) = dStruct.dataTable;
            else
                for gf = 1:numel(loadFile)
                    dStruct = obj.(loadFunc)('FILE',loadFile{gf});
                    if dStruct.error
                        loadStruct.error = 1;
                        return;
                    end
                    % Merge tables
                    if gf == 1
                        obj.data.(loadTable) = dStruct.dataTable;
                    else
                        %by default shift animal ids by 100*(gf-1)
                        dStruct.dataTable.animal_id = dStruct.dataTable.animal_id+100*(gf-1);
                        %                          obj.data.(loadTable) = join( obj.data.(loadTable),dStruct.dataTable,'Keys','animal_id');
                        obj.data.(loadTable) = vertcat( obj.data.(loadTable),dStruct.dataTable);
                    end
                end
            end
            obj.data.(loadTable) = sortrows( obj.data.(loadTable));
        end
        function M = mergeDataFiles(obj,varargin)
            % FUNCTION
            %
            % USAGE:
            % obj.function()
            %
            % ------------------------------------------------------------
            % REQUIRED ARGUMENTS:
            % 'HANDLE_PARENT',parent
            % 'KEY_ID',key_id)
            %
            % OPTIONAL ARGUMENTS:
            %
            % OUTPUTS:
            %
            % NOTES:
            %
            M = struct('error', 0,...
                'dataTable', []);

            validLogical = [0 1];
            p = inputParser;
            p.addParameter('skip',true,@(x) islogical(x) || isnumeric(x));
            p.parse(varargin{:});
            inputs = p.Results;
            %             try
            %                 %check inputs
            %                 p = inputParser;
            %                 p.addParameter('FILE', ...
            %                     [],...
            %                     @(x) ischar(x));
            %                 p.parse(varargin{:});
            %                 inputs = p.Results;
            %             catch err
            %                 obj.status = sprintf('FUNCTION "%s" CALL ERROR: %s', upper(obj.current_function),err.message);
            %                 obj.help('mergeDataFiles');
            %                 return;
            %             end
            obj.status = sprintf('MERGING DATA...');

            % merge viewpoint and grouping -- easy

            % Reconcile group and viewpoint
            % check whether animal ids sets are equivalent
            vp_id = unique(obj.data.viewpoint_data.animal_id);
            g_id = unique(obj.data.group_data.animal_id);

            invp_notgp = setdiff(vp_id,g_id);
            ingp_notvp = setdiff(g_id,vp_id);

            tempGroupDataTable = obj.data.group_data;
            tempViewpointDataTable = obj.data.viewpoint_data;

            if invp_notgp
                if ~inputs.skip
                    obj.status = 'USER INPUT NEEDED!';
                    fprintf('The following animal ids are in the viewpoint data table but not the group table\n')
                    fprintf('%d\n',invp_notgp);
                    qstring = 'Proceed with analysis without these animals?';
                    choice = questdlg(qstring,'Continue analysis',...
                        'Yes','No','No');
                    if strcmp(choice,'No')
                        M.error = 1;
                        obj.status = 'HALTING ANALYSIS';
                        return;
                    end
                end

                % remove animals for merge purposes...
                for i = 1 : numel(invp_notgp)
                    tempViewpointDataTable(tempViewpointDataTable.animal_id == invp_notgp(i),:)= [];
                end
                %                tempViewpointDataTable(tempViewpointDataTable.animal_id == invp_notgp,:)= [];
            end

            if ingp_notvp
                if ~inputs.skip
                    obj.status = 'USER INPUT NEEDED!';
                    fprintf('The following animal ids are in the group data table but not the viewpoint table\n')
                    fprintf('%d\n',ingp_notvp);
                    qstring = 'Proceed with analysis without these animals?';
                    choice = questdlg(qstring,'Continue analysis',...
                        'Yes','No','No');
                    if strcmp(choice,'No')
                        M.error = 1;
                        obj.status = 'HALTING ANALYSIS';
                        return;
                    end
                end

                % remove animals for merge purposes...
                tempGroupDataTable(tempGroupDataTable.animal_id == ingp_notvp,:)= [];
            end

            M.dataTable = join(tempViewpointDataTable,tempGroupDataTable,'Keys','animal_id');

            % merge w/timing -- hard

            % for every row in the merged data table,
            vts = M.dataTable.start;
            eton = obj.data.timing_data.Onset;
            etof = obj.data.timing_data.Offset;


            %             vts=  datetime(M.dataTable.sttime,'Format','HH:mm:SSS');
            %             eton = cellfun(@(x) datetime(x,'Format','HH:mm:SSS'),obj.data.timing_data.onset_date,'UniformOutput',0);
            %             %bump offset down a bit because matlab's "isbetween" is
            %             %inclusive
            %             etof = cellfun(@(x) datetime(x,'Format','HH:mm:SSS')-milliseconds(1),obj.data.timing_data.offset_date,'UniformOutput',0);
            %

            %columns to add to merged table from event timing table
            unNeededCols= {'',...
                'Time',...
                'onset_date',...
                'offset_date',...
                'duration_sec'};
            neededColumns = setdiff(obj.data.timing_data.Properties.VariableNames, unNeededCols);

            % add columns and set to a default value
            for c = 1 :numel(neededColumns)
                % create defaults
                if strcmp(neededColumns{c},'Onset') || strcmp(neededColumns{c},'Offset')
                    M.dataTable.(neededColumns{c}) = repmat(-1,height(M.dataTable),1);
                else

                    M.dataTable.(neededColumns{c}) = repmat({'unspecified'},height(M.dataTable),1);
                end
            end

            %enter timing labels
            %             for e = 1:numel(eton)
            %                 rows = find(isbetween(vts,eton{e},etof{e}));
            %                 for c = 1 :numel(neededColumns)
            %                     M.dataTable.(neededColumns{c})(rows) = obj.data.timing_data.(neededColumns{c})(e);
            %                 end
            %             end
            %
            for e = 1:numel(eton)

                rows = find(vts >= eton(e) & vts < etof(e));
                for c = 1 :numel(neededColumns)
                    try
                        M.dataTable.(neededColumns{c})(rows) = repmat(obj.data.timing_data.(neededColumns{c})(e),numel(rows),1);
                    catch

                    end
                end
            end




            % mark acclimation and common baseline
            neededColumns = {'acclimation' 'common_baseline'};
            %
            %
            %             % add columns and set to a default value
            for c = 1 :numel(neededColumns)
                % create defaults
                M.dataTable.(neededColumns{c}) = repmat([0],height(M.dataTable),1);

            end
            %             t(M.dataTable.xEnd <= obj.meta.acclimationPeriodEnd) = 1
            try
                M.dataTable.acclimation(M.dataTable.xEnd <= obj.meta.acclimationPeriodEnd) = 1;
                %             t((M.dataTable.xEnd <= obj.meta.commonBaselinePeriodEnd) & (M.dataTable.start >= obj.meta.commonBaselinePeriodStart)) = 1;
                M.dataTable.common_baseline((M.dataTable.xEnd <= obj.meta.commonBaselinePeriodEnd) & (M.dataTable.start >= obj.meta.commonBaselinePeriodStart)) = 1;
            catch err
            end
            %             M.dataTable.acclimation(M.dataTable.xEnd <= obj.meta.acclimationPeriodEnd) = 1;

            obj.status = sprintf('MERGING DONE!');

            if isempty(obj.meta.save_location)
                if ismac()
                    obj.meta.save_location = [obj.meta.rootAnalysisOutput '/' obj.meta.project_name];
                else
                    obj.meta.save_location = [obj.meta.rootAnalysisOutput '\' obj.meta.project_name];

                end

            end
            obj.data.merged_data = M.dataTable;

            % %                         mdf = fullfile(obj.meta.save_location,sprintf('mergedData_%s.xlsx',datetime('now','Format','yMMddHHmmSS')));
            % %                         obj.status = sprintf('WRITING MERGE DATA FILE %s...',mdf);
            %                         try
            %                             %                 writetable(M.dataTable,mdf);
            %                             %                 obj.files.merged_file = mdf;
            %                         catch err
            %                             obj.status = sprintf('FUNCTION "%s" MERGE DATA TABLE WRITE ERROR: %s', upper(obj.current_function),err.message);
            %                             M.error = 1;
            %                             return;
            %                         end

        end

        function E = registerNewDataColumn(obj,s,name)
            obj.meta.(s).count = obj.meta.(s).count + 1;
            obj.meta.(s).(name) = obj.meta.(s).count ;

        end
        function G = loadGroupingFile(obj,varargin)
            % FUNCTION
            %
            % USAGE:
            % obj.function()
            %
            % ------------------------------------------------------------
            % REQUIRED ARGUMENTS:
            % 'HANDLE_PARENT',parent
            % 'KEY_ID',key_id)
            %
            % OPTIONAL ARGUMENTS:
            %
            % OUTPUTS:
            %
            % NOTES:
            %


            G = struct('error', 0,...
                'dataTable', []);

            validLogical = [0 1];
            try
                %check inputs
                p = inputParser;
                p.addParameter('FILE', ...
                    [],...
                    @(x) ischar(x));
                p.parse(varargin{:});
                inputs = p.Results;
            catch err
                obj.status = sprintf('FUNCTION "%s" CALL ERROR: %s', upper(obj.current_function),err.message);
                obj.help('loadgroupingFile');
                return;
            end

            obj.status = sprintf('LOADING GROUPING FILE %s...',inputs.FILE);
            if ~exist(   inputs.FILE,'file')
                obj.status = sprintf('ERROR: grouping file could not be loaded: %s',inputs.FILE);
                G.error = 1;
                return;
            end

            obj.status = sprintf('PARSING GROUPING FILE %s...',inputs.FILE);


            T = readtable(inputs.FILE);

            % remove empty trailing column(s)
            idx = false(1,size(T,2));
            for i=1:size(T,2)
                if isnumeric(T{:,i}) && all(isnan(T{:,i}))
                    idx(i) = true;
                end
            end
            T(:,idx) = [];

            %need to figure out alternate header lines
            first_real_row = [];
            for i=1:size(T,2)
               col = find(~isnan(str2double(T{:,i})),1);
               if ~isempty(col) && (isempty(first_real_row) || col < first_real_row)
                  first_real_row = col; 
               end
            end

            altLabels = [];
            if first_real_row ~= 1
                altLabels = strtrim(T{1:first_real_row-1,:});
                T(1:first_real_row-1,:) = [];
            end



            %convert to numbers
            CT = cellfun(@(x) str2num(x),T{:,:},'UniformOutput',0);
            CT = cell2table(CT);
            CT.Properties.VariableNames = T.Properties.VariableNames;

            rows = [];
            for col =  1 : width(CT)
                for row = 1 : height(CT)
                    aid = CT{row,col};
                    if iscell(aid)
                        aid = aid{:};
                    end
                    if ~isempty(aid)

                        row = {aid,CT.Properties.VariableNames{col},col};
                        for al = 1 : size(altLabels,1)
                            row = [row, altLabels{col}];

                        end
                        rows = [rows;row];
                        %                     tempTable = table('animal_id',aid, 'group_label1',CT.Properties.VariableNames{col},'groupNum',col);
                        %                     GT = [GT ; tempTable];
                        %                     if isempty(GT.Properties.VariableNames)
                        %                        pause();
                        %                     end
                    end
                end

            end
            GT = cell2table(rows);
            vn = {'animal_id' 'group_label1' 'groupNum'};
            for al = 1 : size(altLabels,1)
                vn = [vn, sprintf('group_label%d',al+1)];
            end
            GT.Properties.VariableNames = vn;

            % read first line and parse by delimiter -- this gives the the
            % levels of the group


            obj.status = sprintf('GROUPING FILE %s PARSED!',inputs.FILE);

            
            
            
            obj.status = sprintf('RENUMBERING GROUPS!',inputs.FILE);
            
            
            unique_labels  = unique(GT.group_label2) ;
            [bool indx] = ismember(GT.group_label2,unique_labels) ;
            GT.groupNum = indx ;
            
            
            
            obj.status = 'VALIDATING INPUT GROUPING...';
            % do some validation steps

            % make sure the table is well-formatted
            if numel(GT.Properties.VariableNames) ~= 4 || ...
                    ~ismember('animal_id',GT.Properties.VariableNames) || ...
                    ~ismember('group_label1',GT.Properties.VariableNames) || ...
                    ~ismember('group_label2',GT.Properties.VariableNames) || ...
                    ~ismember('groupNum',GT.Properties.VariableNames)
                errorStruct.message = sprintf('Grouping file %s is not formatted properly',inputs.FILE);
                errorStruct.identifier = 'validityCheck:GroupingFile';
                error(errorStruct)
            end
                     
            
            % Some checks to see if group file was sensibly parsed

            if isempty(GT.animal_id)
                obj.status = sprintf('ERROR: No IDs were loaded from file: %s',inputs.FILE);
                G.error = 1;
                return;
            end


            % Make sure that an animal hasn't been assigned to more than
            % 1 group

            u=unique(GT.animal_id);
            n=histc(GT.animal_id,u);
            B = u(n>1);
            if ~isempty(B)
                
                errorStruct.message = [sprintf('The following IDs were found in more than 1 group in the grouping file %s: ',inputs.FILE) sprintf('%d ',B)]
                errorStruct.identifier = 'validityCheck:GroupingFile';
                error(errorStruct)
               
                return;
            end

            % Report any unassigned numbers
            A = 1:numel(GT.animal_id);
            B = setdiff(A,GT.animal_id);

            if ~isempty(B)
                fprintf('NOTE: The following IDs were not found in the grouping file %s:\n',inputs.FILE);
                fprintf('%d\n',B);
                fprintf('\n');
            end
            obj.status = 'INPUTS VALID!';


            G.dataTable = GT;
            obj.status = sprintf('GROUPING FILE %s PARSED',inputs.FILE);
        end
        function T = loadTimingFile(obj,varargin)
            % FUNCTION
            %
            % USAGE:
            % obj.function()
            %
            % ------------------------------------------------------------
            % REQUIRED ARGUMENTS:
            % 'HANDLE_PARENT',parent
            % 'KEY_ID',key_id)
            %
            % OPTIONAL ARGUMENTS:
            %
            % OUTPUTS:
            %
            % NOTES:
            %
            T = struct('error', 0,...
                'dataTable', []);

            validLogical = [0 1];
            try
                %check inputs
                p = inputParser;
                p.addParameter('FILE', ...
                    [],...
                    @(x) ischar(x));
                p.parse(varargin{:});
                inputs = p.Results;
            catch err
                obj.status = sprintf('FUNCTION "%s" CALL ERROR: %s\n\n', upper(obj.current_function),err.message);
                obj.help('loadTimingFile');
                return;
            end

            obj.status = sprintf('LOADING TIMING FILE %s...',inputs.FILE);
            try
                fid = fopen(inputs.FILE,'r');
            catch err
                obj.status = sprintf('ERROR: timing file could not be loaded: %s\n',inputs.FILE);
                T.error = 1;
                return;
            end

            if fid == -1
                obj.status = sprintf('ERROR: timing file could not be loaded: %s\n',inputs.FILE);
                T.error = 1;
                return;
            end

            obj.status = sprintf('PARSING TIMING FILE %s...',inputs.FILE);

            try
                T.dataTable = readtable(inputs.FILE);
            catch
                obj.status = sprintf('ERROR Parsing timing file %s',inputs.FILE);
                T.error = 1;
                return;
            end



            % determine time coding in file
            checkTime = 0;
            needOffset = 0;
            needDuration = 0;
            if ismember('ONSET',upper( T.dataTable.Properties.VariableNames ))
                %onset is recorded
                if ismember('OFFSET',upper( T.dataTable.Properties.VariableNames ))
                    % onset and offset is recorded
                    if ismember('DURATION',upper( T.dataTable.Properties.VariableNames ))
                        % onset, offset, and duration recorded
                        % no additional parameters needed
                    else
                        % duration needed
                        needDuration = 1;
                    end
                else
                    %offset needed
                    if ismember('DURATION',upper( T.dataTable.Properties.VariableNames ))
                        % onset, duration recorded
                        % offset needed
                        needOffset = 1;
                    else
                        % only onset present, not enough time info --
                        % check "time" field
                        checkTime = 1;
                    end
                end
            else
                % onset not recorded, not enough time info --
                % check "time" field
                checkTime = 1;
            end

            if checkTime
                % For now we assume that by default "Time" indicates both the
                % onset and offset of events.  Not an ideal way of representing
                % time

                %first, split the time component based on "-"
                tm = cellfun(@(x) strsplit(x,'-'),T.dataTable.Time,'UniformOutput',0);
                %reshape into # Events x 2 array
                tm = reshape([tm{:}], 2, numel(tm))';
                %remove superfluous whitespace
                tm = strtrim(tm);
                tm = strrep(tm,' ','');

                T.dataTable.onset_date = tm(:,1);
                T.dataTable.offset_date = tm(:,2);
                T.dataTable.duration_sec =  seconds(duration(datetime(tm(:,2)) - datetime(tm(:,1))));

            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Code this once timing is finalized
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if needDuration
            end
            if needOffset
            end







            obj.status = sprintf('TIMING FILE %s PARSED',inputs.FILE);
        end
        function V = loadViewpointFile(obj,varargin)
            % FUNCTION
            %
            % USAGE:
            % obj.function()
            %
            % ------------------------------------------------------------
            % REQUIRED ARGUMENTS:
            % 'HANDLE_PARENT',parent
            % 'KEY_ID',key_id)
            %
            % OPTIONAL ARGUMENTS:
            %
            % OUTPUTS:
            %
            % NOTES:
            %
            V = struct('error', 0,...
                'dataTable', []);

            validLogical = [0 1];
            try
                %check inputs
                p = inputParser;
                p.addParameter('FILE', ...
                    [],...
                    @(x) ischar(x));
                p.parse(varargin{:});
                inputs = p.Results;
            catch err
                obj.status = sprintf('FUNCTION "%s" CALL ERROR: %s', upper(obj.current_function),err.message);
                obj.help('loadViewpointFile');
                return;
            end

            obj.status = sprintf('LOADING VIEWPOINT FILE %s...',inputs.FILE);
            %             try
            %                 fid = fopen(inputs.FILE,'r');
            %             catch err
            %                 obj.status = sprintf('ERROR: Viewpoint file could not be loaded: %s',inputs.FILE);
            %                 V.error = 1;
            %                 return;
            %             end
            %
            %             if fid == -1
            %                 obj.status = sprintf('ERROR: Viewpoint file could not be loaded: %s',inputs.FILE);
            %                 V.error = 1;
            %                 return;
            %             end
            obj.status = 'PARSING MAY TAKE SOME TIME...';
            obj.status = sprintf('PARSING VIEWPOINT FILE %s...',inputs.FILE);

            
            [~,~,fext] = fileparts(inputs.FILE);
            switch fext
                case {'.txt','.XLS'}
                    V.dataTable = readtable(inputs.FILE,'FileType','text','Delimiter','\t');
                otherwise
                    V.dataTable = readtable(inputs.FILE);
            end
            
            if size(V.dataTable,1) == 0
                errorStruct.message = sprintf('Viewpoint file %s possibly corrupt, no data read in',inputs.FILE) ;
                errorStruct.identifier = 'validityCheck:viewpoint';
                error(errorStruct)
            end
            
            V.dataTable.sttime = cellstr(datestr(V.dataTable.sttime,'HH:MM:SS')); % match .txt input, otherwise combine error later
            V.dataTable.stdate = cellstr(datestr(V.dataTable.stdate,'dd/mm/yyyy')); % match .txt input, otherwise combine error later



            obj.status = 'MAPPING ANIMAL LOCATIONS TO NUMERICAL...';

            V.dataTable.animal_id = cellfun(@(x) str2num(strrep(x,'w','')),V.dataTable.location);

            % obj.status = 'REFORMATTING TIMESTAMP REPRESENTATIONS...';
            % idiosyncratic format manipulations
            %V.dataTable.sttime = cellfun(@(x) datestr(x,'HH:MM:SS'),num2cell(V.dataTable.sttime),'UniformOutput',0);

            obj.status = sprintf('VIEWPOINT FILE %s PARSED',inputs.FILE);
        end
        function S = computeDerivedVariables2(obj)

            %             1.	STIM_mean: 1 way
            %             2.	POST_mean: 1 way
            %             3.	S1_mean: 1 way
            %             4.	S2_mean: 1 way
            %             5.	S3_mean: 1 way
            %             6.	S4_mean: 1 way
            %             7.	S5_mean: 1 way
            %             8.	P1_std: 1 way
            %             9.	P2_std: 1 way
            %             10.	P3_std: 1 way
            %             11.	P4_std: 1 way
            %             12.	P5_std: 1 way



            obj.statistics.uni.STIM_mean_nrm = obj.statistics.uni.STIM_mean - obj.statistics.uni.BASELINE_mean ;
            obj.statistics.uni.POST_mean_nrm = obj.statistics.uni.POST_mean - obj.statistics.uni.BASELINE_mean ;
            obj.statistics.uni.POST_std_nrm = obj.statistics.uni.POST_std - obj.statistics.uni.BASELINE_std ;

            obj.statistics.uni.S1_mean_nrm = obj.statistics.uni.S1_mean - obj.statistics.uni.BASELINE_mean ;
            obj.statistics.uni.S2_mean_nrm = obj.statistics.uni.S2_mean - obj.statistics.uni.BASELINE_mean ;
            obj.statistics.uni.S3_mean_nrm = obj.statistics.uni.S3_mean - obj.statistics.uni.BASELINE_mean ;
            obj.statistics.uni.S4_mean_nrm = obj.statistics.uni.S4_mean - obj.statistics.uni.BASELINE_mean ;
            obj.statistics.uni.S5_mean_nrm = obj.statistics.uni.S5_mean - obj.statistics.uni.BASELINE_mean ;

            obj.statistics.uni.P1_mean_nrm = obj.statistics.uni.P1_mean - obj.statistics.uni.BASELINE_mean ;
            obj.statistics.uni.P2_mean_nrm = obj.statistics.uni.P2_mean - obj.statistics.uni.BASELINE_mean ;
            obj.statistics.uni.P3_mean_nrm = obj.statistics.uni.P3_mean - obj.statistics.uni.BASELINE_mean ;
            obj.statistics.uni.P4_mean_nrm = obj.statistics.uni.P4_mean - obj.statistics.uni.BASELINE_mean ;
            obj.statistics.uni.P5_mean_nrm = obj.statistics.uni.P5_mean - obj.statistics.uni.BASELINE_mean;



            %             13.	STIM-POST : 2 way, derived
            %             14.	S1-S5_mean : 2 way, derived
            %             15.	P1-P5_mean : 2 way, derived
            obj.statistics.uni.BASELINEvsSTIM_mean = obj.statistics.uni.STIM_mean - obj.statistics.uni.BASELINE_mean ;
            obj.statistics.uni.STIMvPOST_mean = obj.statistics.uni.STIM_mean - obj.statistics.uni.POST_mean ;
            obj.statistics.uni.S1vS5_mean = obj.statistics.uni.S1_mean - obj.statistics.uni.S5_mean ;

        end
        function univariateStats2(obj)



            %split stats by group
            groups = unique(obj.statistics.uni.group_label2) ;

            NumGroups = numel(groups) ;


            for g = 1 :  NumGroups
                rows = ismember(obj.statistics.uni.group_label2,groups{g}) ;
                uni_group{g} = obj.statistics.uni(rows,:);
            end

            %             One-way stats
            obj.onewayANOVA('VAR','BASELINE_mean','GRP_COL','group_label2', 'TITLE', 'BASELINE_mean') ;
            obj.onewayANOVA('VAR','BASELINE_conmean','GRP_COL','group_label2', 'TITLE', 'BASELINE_conmean') ;
            obj.onewayANOVA('VAR','BASELINE_std','GRP_COL','group_label2', 'TITLE', 'BASELINE_std') ;
            obj.onewayANOVA('VAR','STIM_mean','GRP_COL','group_label2', 'TITLE', 'STIM_mean') ;
            obj.onewayANOVA('VAR','STIM_conmean','GRP_COL','group_label2', 'TITLE', 'STIM_conmean') ;
            obj.onewayANOVA('VAR','STIM_std','GRP_COL','group_label2', 'TITLE', 'STIM_std') ;
            obj.onewayANOVA('VAR','POST_mean','GRP_COL','group_label2', 'TITLE', 'POST_mean') ;
            obj.onewayANOVA('VAR','POST_conmean','GRP_COL','group_label2', 'TITLE', 'POST_conmean') ;
            obj.onewayANOVA('VAR','POST_std','GRP_COL','group_label2', 'TITLE', 'POST_std') ;
            obj.onewayANOVA('VAR','S1_mean','GRP_COL','group_label2', 'TITLE', 'S1_mean') ;
            obj.onewayANOVA('VAR','S1_conmean','GRP_COL','group_label2', 'TITLE', 'S1_conmean') ;
            obj.onewayANOVA('VAR','P1_mean','GRP_COL','group_label2', 'TITLE', 'P1_mean') ;
            obj.onewayANOVA('VAR','P1_std','GRP_COL','group_label2', 'TITLE', 'P1_std') ;
            obj.onewayANOVA('VAR','S2_mean','GRP_COL','group_label2', 'TITLE', 'S2_mean') ;
            obj.onewayANOVA('VAR','S2_conmean','GRP_COL','group_label2', 'TITLE', 'S2_conmean') ;
            obj.onewayANOVA('VAR','P2_mean','GRP_COL','group_label2', 'TITLE', 'P2_mean') ;
            obj.onewayANOVA('VAR','P2_std','GRP_COL','group_label2', 'TITLE', 'P2_std') ;
            obj.onewayANOVA('VAR','S3_mean','GRP_COL','group_label2', 'TITLE', 'S3_mean') ;
            obj.onewayANOVA('VAR','S3_conmean','GRP_COL','group_label2', 'TITLE', 'S3_conmean') ;
            obj.onewayANOVA('VAR','P3_mean','GRP_COL','group_label2', 'TITLE', 'P3_mean') ;
            obj.onewayANOVA('VAR','P3_std','GRP_COL','group_label2', 'TITLE', 'P3_std') ;
            obj.onewayANOVA('VAR','S4_mean','GRP_COL','group_label2', 'TITLE', 'S4_mean') ;
            obj.onewayANOVA('VAR','S4_conmean','GRP_COL','group_label2', 'TITLE', 'S4_conmean') ;
            obj.onewayANOVA('VAR','P4_mean','GRP_COL','group_label2', 'TITLE', 'P4_mean') ;
            obj.onewayANOVA('VAR','P4_std','GRP_COL','group_label2', 'TITLE', 'P4_std') ;
            obj.onewayANOVA('VAR','S5_mean','GRP_COL','group_label2', 'TITLE', 'S5_mean') ;
            obj.onewayANOVA('VAR','S5_conmean','GRP_COL','group_label2', 'TITLE', 'S5_conmean') ;
            obj.onewayANOVA('VAR','P5_mean','GRP_COL','group_label2', 'TITLE', 'P5_mean') ;
            obj.onewayANOVA('VAR','P5_std','GRP_COL','group_label2', 'TITLE', 'P5_std') ;



            obj.onewayCHI('var','BASELINE_propMv','GRP_COL','group_label2', 'TITLE', 'BASELINE_propMv') ;
            obj.onewayBAR('var','STIM_propMv','GRP_COL','group_label2', 'TITLE', 'STIM_propMv') ;
            obj.onewayBAR('var','POST_propMv','GRP_COL','group_label2', 'TITLE', 'POST_propMv') ;
            obj.onewayBAR('var','S1_propMv','GRP_COL','group_label2', 'TITLE', 'S1_propMv') ;
            obj.onewayBAR('var','P1_propMv','GRP_COL','group_label2', 'TITLE', 'P1_propMv') ;
            obj.onewayBAR('var','S2_propMv','GRP_COL','group_label2', 'TITLE', 'S2_propMv') ;
            obj.onewayBAR('var','P2_propMv','GRP_COL','group_label2', 'TITLE', 'P2_propMv') ;
            obj.onewayBAR('var','S3_propMv','GRP_COL','group_label2', 'TITLE', 'S3_propMv') ;
            obj.onewayBAR('var','P3_propMv','GRP_COL','group_label2', 'TITLE', 'P3_propMv') ;
            obj.onewayBAR('var','S4_propMv','GRP_COL','group_label2', 'TITLE', 'S4_propMv') ;
            obj.onewayBAR('var','P4_propMv','GRP_COL','group_label2', 'TITLE', 'P4_propMv') ;
            obj.onewayBAR('var','S5_propMv','GRP_COL','group_label2', 'TITLE', 'S5_propMv') ;
            obj.onewayBAR('var','P5_propMv','GRP_COL','group_label2', 'TITLE', 'P5_propMv') ;


            obj.onewayBAR('var','S1_premove','GRP_COL','group_label2', 'TITLE', 'S1_premove') ;
            obj.onewayBAR('var','S2_premove','GRP_COL','group_label2', 'TITLE', 'S2_premove') ;
            obj.onewayBAR('var','S3_premove','GRP_COL','group_label2', 'TITLE', 'S3_premove') ;
            obj.onewayBAR('var','S4_premove','GRP_COL','group_label2', 'TITLE', 'S4_premove') ;
            obj.onewayBAR('var','S5_premove','GRP_COL','group_label2', 'TITLE', 'S5_premove') ;

            obj.twowayANOVA('VAR1','BASELINE_mean', 'VAR2','STIM_mean','GRP_COL','group_label2','TITLE', 'base_vs_STIM') ;
            obj.twowayANOVA('VAR1','BASELINE_mean', 'VAR2','POST_mean','GRP_COL','group_label2','TITLE', 'base_vs_POST') ;
            obj.twowayANOVA('VAR1','BASELINE_std', 'VAR2','POST_std','GRP_COL','group_label2','TITLE', 'base_vs_post_std') ;
            obj.twowayANOVA('VAR1','BASELINE_mean', 'VAR2','S1_mean','GRP_COL','group_label2','TITLE', 'base_vs_S1') ;
            obj.twowayANOVA('VAR1','BASELINE_mean', 'VAR2','S2_mean','GRP_COL','group_label2','TITLE', 'base_vs_S2') ;
            obj.twowayANOVA('VAR1','BASELINE_mean', 'VAR2','S3_mean','GRP_COL','group_label2','TITLE', 'base_vs_S3') ;
            obj.twowayANOVA('VAR1','BASELINE_mean', 'VAR2','S4_mean','GRP_COL','group_label2','TITLE', 'base_vs_S4') ;
            obj.twowayANOVA('VAR1','BASELINE_mean', 'VAR2','S5_mean','GRP_COL','group_label2','TITLE', 'base_vs_S5') ;
            obj.twowayANOVA('VAR1','BASELINE_mean', 'VAR2','P1_mean','GRP_COL','group_label2','TITLE', 'base_vs_P1') ;
            obj.twowayANOVA('VAR1','BASELINE_mean', 'VAR2','P2_mean','GRP_COL','group_label2','TITLE', 'base_vs_P2') ;
            obj.twowayANOVA('VAR1','BASELINE_mean', 'VAR2','P3_mean','GRP_COL','group_label2','TITLE', 'base_vs_P3') ;
            obj.twowayANOVA('VAR1','BASELINE_mean', 'VAR2','P4_mean','GRP_COL','group_label2','TITLE', 'base_vs_P4') ;
            obj.twowayANOVA('VAR1','BASELINE_mean', 'VAR2','P5_mean','GRP_COL','group_label2','TITLE', 'base_vs_P5') ;
            obj.twowayANOVA('VAR1','STIM_mean', 'VAR2','POST_mean','GRP_COL','group_label2','TITLE', 'STIM_vs_POST') ;
            obj.twowayANOVA('VAR1','S1_mean', 'VAR2','S5_mean','GRP_COL','group_label2','TITLE', 'S1_vs_S5') ;




        end
        function stats = computeAnimalWiseProperties_arch(obj, varargin)
            %sort the merged data table by time ascending
            obj.data.merged_data = sortrows(obj.data.merged_data, {'animal_id' 'start'}) ;
            % attach unique row number
            obj.data.merged_data.rowNum = [1 : height(obj.data.merged_data)]' ;

            %set up row sets  for later searches
            BASELINE_ROWS = obj.data.merged_data.rowNum(obj.data.merged_data.common_baseline==1);
            STIM_ROWS = obj.data.merged_data.rowNum(obj.data.merged_data.event_period==1);
            POST_ROWS = obj.data.merged_data.rowNum(obj.data.merged_data.event_postperiod==1);

            for events = 1 : obj.meta.NumEvent
                EV{events} = obj.data.merged_data.rowNum(obj.data.merged_data.event_iteration==events);
            end
            all_animals = unique(obj.data.merged_data.animal_id) ;

            %set up stat structure

            stats.uni = table(all_animals,'VariableNames', {'animal_id'}) ;
            stats.ts = table(all_animals,'VariableNames', {'animal_id'}) ;
            stats.uni = join(stats.uni,obj.data.group_data,'Keys','animal_id') ;
            stats.ts = join(stats.ts,obj.data.group_data,'Keys','animal_id') ;

            %             stats.uni.BASELINE_mean = repmat(0,height(stats.uni),1);
            %             stats.uni.BASELINE_std = repmat(0,height(stats.uni),1);


            % for each animal compute:
            for a = 1 : numel(all_animals)
                thisAnimal = all_animals(a) ;
                thisRow = find(stats.uni.animal_id == thisAnimal);
                an = obj.data.merged_data.rowNum(obj.data.merged_data.animal_id==thisAnimal);

                br = intersect(an, BASELINE_ROWS) ;
                s = obj.computeBasicStats('ROWS',br,'DEP_VAR',obj.meta.dep_var) ;
                % baseline
                stats.uni.BASELINE_mean(thisRow) = s.mean ;
                stats.uni.BASELINE_std(thisRow) = s.std ;
                stats.uni.BASELINE_propMv(thisRow) = s.propMv ;
                stats.uni.BASELINE_conmean(thisRow) = s.conmean ;
                stats.uni.BASELINE_constd(thisRow) = s.constd ;
                stats.ts.BASELINE_ERA(thisRow) = {s.ERA'} ;
                stats.ts.BASELINE_bERA(thisRow) = {s.bERA'} ;
                stats.ts.BASELINE_cERA(thisRow) = {s.cERA'} ;
                %                 stats.ts.BASELINE_ERASM(thisRow) = {s.ERA_SM'} ;


                %                 % for all events + postevents together
                %                 br = intersect(an, union(STIM_ROWS,POST_ROWS)) ;
                % %                 s = obj.computeBasicStats('ROWS',br,'DEP_VAR',obj.meta.dep_var,'BY_EVENT',1) ;
                % %                 stats.uni.allepe_mean(thisRow) = s.mean ;
                % %                 stats.uni.allepe_std(thisRow) = s.std ;
                %                  stats.ts.STIMPOST_ERA(thisRow) = {s.ERA'} ;
                % %                  stats.ts.STIMPOST_ERASM(thisRow) = {s.ERA_SM'} ;

                % STIM
                br = intersect(an, STIM_ROWS) ;
                s = obj.computeBasicStats('ROWS',br,'DEP_VAR',obj.meta.dep_var,'BY_EVENT',0) ;
                stats.uni.STIM_mean(thisRow) = s.mean ;
                stats.uni.STIM_propMv(thisRow) = s.propMv ;
                stats.uni.STIM_std(thisRow) = s.std ;
                stats.uni.STIM_conmean(thisRow) = s.conmean ;
                stats.uni.STIM_constd(thisRow) = s.constd ;


                %                 stats.uni.STIM_std(thisRow) = s.std ;
                %                 stats.ts.STIM_ERA(thisRow) = {s.ERA'} ;
                %                 stats.ts.STIM_ERASM(thisRow) = {s.ERA_SM'} ;

                % POST
                br = intersect(an, POST_ROWS) ;
                s = obj.computeBasicStats('ROWS',br,'DEP_VAR',obj.meta.dep_var,'BY_EVENT',0) ;
                stats.uni.POST_mean(thisRow) = s.mean ;
                stats.uni.POST_propMv(thisRow) = s.propMv ;
                stats.uni.POST_std(thisRow) = s.std ;
                stats.uni.POST_conmean(thisRow) = s.conmean ;
                stats.uni.POST_constd(thisRow) = s.constd ;

                %                 stats.ts.POST_ERA(thisRow) = {s.ERA'} ;
                %                 stats.ts.POST_ERASM(thisRow) = {s.ERA_SM'} ;

                %                 %                 STIM + POST (for ERA only
                %                 br = intersect(an, union(STIM_ROWS,POST_ROWS)) ;
                %                 s = obj.computeBasicStats('ROWS',br,'DEP_VAR',obj.meta.dep_var,'BY_EVENT',1) ;
                %                 %                                 stats.uni.allepe_mean(thisRow) = s.mean ;
                %                 %                                 stats.uni.allepe_std(thisRow) = s.std ;
                %                 stats.ts.STIMPOST_ERA(thisRow) = {s.ERA'} ;
                % %                  stats.ts.STIMPOST_ERASM(thisRow) = {s.ERA_SM'} ;

                for e = 1 : obj.meta.NumEvent
                    br = intersect(an, intersect(STIM_ROWS,EV{e})) ;
                    s = obj.computeBasicStats('ROWS',br,'DEP_VAR',obj.meta.dep_var,'BY_EVENT',0) ;
                    stats.uni.(genvarname(sprintf('S%d_mean',e)))(thisRow) = s.mean ;
                    stats.uni.(genvarname(sprintf('S%d_std',e)))(thisRow) = s.std ;
                    stats.uni.(genvarname(sprintf('S%d_propMv',e)))(thisRow) = s.propMv ;
                    stats.uni.(genvarname(sprintf('S%d_premove',e)))(thisRow) = s.premove ;
                    stats.uni.(genvarname(sprintf('S%d_conmean',e)))(thisRow) = s.conmean ;
                    stats.uni.(genvarname(sprintf('S%d_constd',e)))(thisRow) = s.constd ;



                    br = intersect(an, intersect(POST_ROWS,EV{e})) ;
                    s = obj.computeBasicStats('ROWS',br,'DEP_VAR',obj.meta.dep_var,'BY_EVENT',0) ;
                    stats.uni.(genvarname(sprintf('P%d_mean',e)))(thisRow) = s.mean ;
                    stats.uni.(genvarname(sprintf('P%d_std',e)))(thisRow) = s.std ;
                    stats.uni.(genvarname(sprintf('P%d_propMv',e)))(thisRow) = s.propMv ;
                    stats.uni.(genvarname(sprintf('P%d_premove',e)))(thisRow) = s.premove ;
                    stats.uni.(genvarname(sprintf('P%d_conmean',e)))(thisRow) = s.conmean ;
                    stats.uni.(genvarname(sprintf('P%d_constd',e)))(thisRow) = s.constd ;

                    br = intersect(an, intersect(union(STIM_ROWS,POST_ROWS),EV{e})) ;
                    s = obj.computeBasicStats('ROWS',br,'DEP_VAR',obj.meta.dep_var,'BY_EVENT',0) ;
                    stats.ts.(genvarname(sprintf('STIMPOST_S%d_ERA',e)))(thisRow) = {s.ERA'} ;
                    stats.ts.(genvarname(sprintf('STIMPOST_S%d_bERA',e)))(thisRow) = {s.bERA'} ;
                    stats.ts.(genvarname(sprintf('STIMPOST_S%d_cERA',e)))(thisRow) = {s.cERA'} ;
                    %
                end
            end
            obj.statistics = stats;

            %             1.	STIM_mean: 1 way
            %             2.	POST_mean: 1 way
            %             3.	S1_mean: 1 way
            %             4.	S2_mean: 1 way
            %             5.	S3_mean: 1 way
            %             6.	S4_mean: 1 way
            %             7.	S5_mean: 1 way
            %             8.	P1_std: 1 way
            %             9.	P2_std: 1 way
            %             10.	P3_std: 1 way
            %             11.	P4_std: 1 way
            %             12.	P5_std: 1 way
            %             13.	STIM-POST : 2 way, derived
            %             14.	S1-S5_mean : 2 way, derived
            %             15.	P1-P5_mean : 2 way, derived

            obj.computeDerivedVariables() ;
            obj.computeTSparams() ;

            obj.univariateStats2 ;
            obj.printStatTable ;


        end

    end
    methods (Static = true, Access = protected)

        function    caller_name = current_function( inarg )
            %
            %   See also: mfilename
            dbk = dbstack( 1 );
            if isempty( dbk )
                str = 'base';
            else
                str = dbk(1).name;
            end
            ixf = find( str == '.', 1, 'first' );
            if isempty( ixf ) || ( nargin==1 && strcmp( inarg, '-full' ) )
                caller_name = str;
            else
                caller_name = str( ixf+1 : end );
            end
        end
        function    initArt(v)

            fprintf('   ,-----.\n  /____  /\n /     \\/ _\n/ O     \\/ |\n>          |\n\\ <)    /\\_|\n \\_____/\\    \n  \\      \\   \n   `-----''\n\n');

            mfile = mfilename('fullpath');
            mfile = fileparts(mfile);
            %             sl = strfind(mfile,filesep);
            %             mfile = mfile(1:sl(end-1));

            %             obj.meta_code_locations = {fullfile(mfile, '\parentClass'),...
            %                 fullfile(mfile,'\subclasses'),...
            %                 fullfile(mfile,'\subclasses\individualBehavioralClasses')};

            % count lines of code  -- mostly for fun


            fprintf('ANIMAL BEHAVIOR PARSER VERSION %s\n',v);
            fprintf('8/22/2016 -Brent Vander Wyk\n');


        end
        function printTable(fid,T,cellCols)
            fprintf(fid,'\n');
            % print variable names
            fprintf(fid,'%s\t',T.Properties.VariableNames{:});
            fprintf(fid,'\n');

            for i = 1 : height(T)
                if cellCols ==1
                    fprintf(fid,'%s\t%d\t%.2f\n',T{i,1}{:},T{i,2},T{i,3}) ;
                elseif cellCols ==2
                    fprintf(fid,'%s\t%s\t%d\t%.2f\n',T{i,1}{:},T{i,2}{:},T{i,3},T{i,4}) ;
                else
                    fprintf(fid,'%s\t%s\t%d\t%d\t%.2f\n',T{i,1}{:},T{i,2}{:},T{i,3},T{i,4},T{i,5}) ;
                end


            end
            fprintf(fid,'\n');


        end
        function data = computeDVSummary(data)

            if isfield(data,'group_count')
                groups = fieldnames(data.data);

                % check if there are post-periods or just events (assuming
                % all groups have the same event types)

                %                 levels = unique(data.data.(groups{1}).events);
                levels = unique(data.data.(groups{1}).events(~cellfun('isempty',data.data.(groups{1}).events)));

                data.data.event_levels = levels;
                data.data.group_means_byevent = zeros(numel(levels), data.group_count);
                data.data.group_std_byevent = zeros(numel(levels), data.group_count);
                data.data.group_sem_byevent = zeros(numel(levels), data.group_count);
                data.data.group_cnts = zeros(1, data.group_count);



                for g =  1 : data.group_count

                    data.data.group_cnts(g)=  data.data.(groups{g}).cnt;
                    data.data.(groups{g}).means_by_animal = zeros(numel(levels),data.data.(groups{g}).cnt);
                    for thisLevel = 1 : numel(levels)
                        theseRows = find(ismember(data.data.(groups{g}).events(:,1),levels(thisLevel)));

                        data.data.(groups{g}).means_by_animal(thisLevel,:) = mean(data.data.(groups{g}).DV(theseRows,:),1);
                        data.data.group_means_byevent(thisLevel,g) = mean(data.data.(groups{g}).means_by_animal(thisLevel,:));
                        data.data.group_std_byevent(thisLevel,g) = std(data.data.(groups{g}).means_by_animal(thisLevel,:));

                        data.data.group_sem_byevent(thisLevel,g) = data.data.group_std_byevent(thisLevel,g)/sqrt(data.data.group_cnts(g));
                    end
                end

            else
                %no group

            end
        end
        function data = datawiseSubtraction(d1,d2,groups)
            if groups
                groupList = fieldnames(d1.data);
                for g =  1 : d1.group_count
                    thisGroup = groupList{g};
                    try
                        d1.data.(thisGroup).DV =  d1.data.(thisGroup).DV - repmat(d2.data.(thisGroup).means_by_animal,size(d1.data.(thisGroup).DV,1),1);
                    catch err
                    end
                end


            else
                % only handles groups right now
            end


            data = d1;

        end
    end


    methods (Static)

        function t = timeStitch(dat,a)
            try
                t = horzcat(dat.BASELINE_ERA{a}, dat.STIMPOST_S1_ERA{a},dat.STIMPOST_S2_ERA{a},dat.STIMPOST_S3_ERA{a},dat.STIMPOST_S4_ERA{a},dat.STIMPOST_S5_ERA{a});
            catch err
                t = [];
            end
        end
        function T = readTimingFile(varargin)
            % FUNCTION
            %
            % USAGE:
            % obj.function()
            %
            % ------------------------------------------------------------
            % REQUIRED ARGUMENTS:
            % 'HANDLE_PARENT',parent
            % 'KEY_ID',key_id)
            %
            % OPTIONAL ARGUMENTS:
            %
            % OUTPUTS:
            %
            % NOTES:
            %
            T = struct('error', 0,...
                'dataTable', []);

            try
                %check inputs
                p = inputParser;
                p.addParameter('FILE', ...
                    [],...
                    @(x) ischar(x));
                p.parse(varargin{:});
                inputs = p.Results;
            catch err
                fprintf('READ_TIMING_FILE ERROR: %s\n\n',err.message);
                return
            end

            obj.status = sprintf('LOADING TIMING FILE %s...',inputs.FILE);
            try
                fid = fopen(inputs.FILE,'r');
            catch
                fprintf('ERROR: timing file could not be loaded: %s\n',inputs.FILE);
                T.error = 1;
                return
            end

            if fid == -1
                fprintf('ERROR: timing file could not be loaded: %s\n',inputs.FILE);
                T.error = 1;
                return
            end

            %fprintf('PARSING TIMING FILE %s...',inputs.FILE);

            try
                T.dataTable = readtable(inputs.FILE);
            catch
                fprintf('ERROR Parsing timing file %s',inputs.FILE);
                T.error = 1;
                return
            end



            % determine time coding in file
            checkTime = 0;
            needOffset = 0;
            needDuration = 0;
            if ismember('ONSET',upper( T.dataTable.Properties.VariableNames ))
                %onset is recorded
                if ismember('OFFSET',upper( T.dataTable.Properties.VariableNames ))
                    % onset and offset is recorded
                    if ismember('DURATION',upper( T.dataTable.Properties.VariableNames ))
                        % onset, offset, and duration recorded
                        % no additional parameters needed
                    else
                        % duration needed
                        needDuration = 1;
                    end
                else
                    %offset needed
                    if ismember('DURATION',upper( T.dataTable.Properties.VariableNames ))
                        % onset, duration recorded
                        % offset needed
                        needOffset = 1;
                    else
                        % only onset present, not enough time info --
                        % check "time" field
                        checkTime = 1;
                    end
                end
            else
                % onset not recorded, not enough time info --
                % check "time" field
                checkTime = 1;
            end

            if checkTime
                % For now we assume that by default "Time" indicates both the
                % onset and offset of events.  Not an ideal way of representing
                % time

                %first, split the time component based on "-"
                tm = cellfun(@(x) strsplit(x,'-'),T.dataTable.Time,'UniformOutput',0);
                %reshape into # Events x 2 array
                tm = reshape([tm{:}], 2, numel(tm))';
                %remove superfluous whitespace
                tm = strtrim(tm);
                tm = strrep(tm,' ','');

                T.dataTable.onset_date = tm(:,1);
                T.dataTable.offset_date = tm(:,2);
                T.dataTable.duration_sec =  seconds(duration(datetime(tm(:,2)) - datetime(tm(:,1))));

            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Code this once timing is finalized
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if needDuration
            end
            if needOffset
            end


            fclose(fid) ;




            %fprintf('TIMING FILE %s PARSED',inputs.FILE);
        end
        function G = readGroupingFile(FILE,varargin)
            % Function ripped from protected method loadGroupingFile for
            % the purpose of making this code available to external
            % scripts. (loading these files is annoying)

            % Parse inputs
            p = inputParser;
            p.addRequired('FILE',@(x) ischar(x) && exist(x,'file'));
            p.parse(FILE,varargin{:});
            inputs = p.Results;

            T = readtable(inputs.FILE);

            %need to figure out alternate header lines

            col1 = T{:,1};
            convcol1 = cellfun(@(x) str2num(x),col1,'UniformOutput',0);
            first_real_row = find(cellfun(@(x) ~isempty(x),convcol1),1);

            altLabels = [];
            if first_real_row ~= 1
                altLabels = strtrim(T{1:first_real_row-1,:});
                T(1:first_real_row-1,:) = [];
            end

            %convert to numbers
            CT = cellfun(@(x) str2num(x),T{:,:},'UniformOutput',0);
            CT = cell2table(CT);
            CT.Properties.VariableNames = T.Properties.VariableNames;

            rows = [];
            for col =  1 : width(CT)
                for row = 1 : height(CT)
                    aid = CT{row,col};
                    if iscell(aid)
                        aid = aid{:};
                    end
                    if ~isempty(aid)

                        row = {aid,CT.Properties.VariableNames{col},col};
                        for al = 1 : size(altLabels,1)
                            row = [row, altLabels{col}];

                        end
                        rows = [rows;row];
                        %                     tempTable = table('animal_id',aid, 'group_label1',CT.Properties.VariableNames{col},'groupNum',col);
                        %                     GT = [GT ; tempTable];
                        %                     if isempty(GT.Properties.VariableNames)
                        %                        pause();
                        %                     end
                    end
                end

            end
            GT = cell2table(rows);
            vn = {'animal_id' 'group_label1' 'groupNum'};
            for al = 1 : size(altLabels,1)
                vn = [vn, sprintf('group_label%d',al+1)];
            end
            GT.Properties.VariableNames = vn;

            % read first line and parse by delimiter -- this gives the the
            % levels of the group

            % do some validation steps

            % Some checks to see if group file was sensibly parsed

            if isempty(GT.animal_id)
                fprintf('ERROR: No IDs were loaded from file: %s\n',inputs.FILE);
                G.error = 1;
                return;
            end


            % Make sure that an animal hasn't been assigned to more than
            % 1 group

            u=unique(GT.animal_id);
            n=histc(GT.animal_id,u);
            B = u(n>1);
            if ~isempty(B)

                fprintf('The following IDs were found in more than 1 group in the grouping file %s\n',inputs.FILE);
                fprintf('%d\n',B);
                fprintf('\n');
                G.error = 1;
                return
            end

            % Report any unassigned numbers
            A = 1:numel(GT.animal_id);
            B = setdiff(A,GT.animal_id);

            if ~isempty(B)
                fprintf('NOTE: The following IDs were not found in the grouping file %s:\n',inputs.FILE);
                fprintf('%d\n',B);
                fprintf('\n');
            end
            %obj.status = 'INPUTS VALID!';
            % Create a table
            % [~,columnHeader,~] = fileparts(inputs.FILE);
            %             if t > 1
            %                 columnHeader = [columnHeader sprintf('_alt%d',t-1)];
            %             end

            % Loop over levels

            %             for L = 1 : size(groupLevels,2)
            %
            %
            %                 newTable = table(assignmentArray{L}(:));
            %                 newTable.Properties.VariableNames = {'animal_id'};
            %                 % Loop over names and add to table
            %                 for N = 1 : size(groupLevels,1)
            %
            %
            %                     % Main group variable name
            % %                     [~,columnHeader,~] = fileparts(inputs.FILE);
            %                     columnHeader = 'GroupLabel_1';
            %                     % Alternate group level names
            %                     if N > 1
            % %                         columnHeader = [columnHeader sprintf('_alt%d',N-1)];
            % %                             columnHeader = sprintf('Grouplabel_%d',N);
            %                     end
            %                     columnHeader = ['grp_' columnHeader];
            %                     newTable.(columnHeader) = repmat(groupLevels(N,L),height(newTable),1);
            %                 end
            %                 if L == 1
            %                     fullTable = newTable;
            %                 else
            %                     fullTable = vertcat(fullTable,newTable);
            %                 end
            %                     %%
            %             end

            G.dataTable = GT;
            %obj.status = sprintf('GROUPING FILE %s PARSED',inputs.FILE);
        end
    end
end
