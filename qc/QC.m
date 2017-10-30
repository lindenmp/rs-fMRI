% ------------------------------------------------------------------------------
% This script performs the primary analyses for
% Parkes, Fulcher, Yucel, Fornito. An evaluation of the efficacy, reliability, 
% and sensitivity of motion correction strategies for resting-state functional MRI
% 
% Copyright (C) 2017, Linden Parkes <lindenparkes@gmail.com>,
% ------------------------------------------------------------------------------
clear all; close all; clc

% ------------------------------------------------------------------------------
% Set string switches
% ------------------------------------------------------------------------------
Projects = {'Beijing_Zang','UCLA','OCDPG','NYU_2','COBRE','COBRE_HCTSA','UCLA_HCTSA','NAMIC_HCTSA','ML_SNS'};
WhichProject = Projects{2};

WhichParc = 'Gordon'; % 'Gordon' 'Power'

if ismember('OCDPG',WhichProject,'rows') | ismember('UCLA',WhichProject,'rows') | ismember('COBRE',WhichProject,'rows')
	WhichSplit = 'Diagnostic'; % 'Motion' 'Diagnostic'
	% WhichSplit = 'Motion'; % 'Motion' 'Diagnostic'
	% Note, this only effect the NBS parts of the script.
	% 'Group' variable ALWAYS represents case-control (e.g., 1 = HC, 2 = patients)
elseif ismember('Beijing_Zang',WhichProject,'rows')
	WhichSplit = 'Motion'; % Only motion for Beijing
end

% ------------------------------------------------------------------------------
% Set logical switches
% ------------------------------------------------------------------------------
runPlot = false;
runOverlapPlots = false;

fprintf(1, 'Running rfMRI QC.\n\tDataset: %s\n\tParcelation: %s\n',WhichProject,WhichParc);
if ismember('OCDPG',WhichProject,'rows') | ismember('UCLA',WhichProject,'rows') | ismember('COBRE',WhichProject,'rows')
	fprintf(1,'\tNBS: %s\n',WhichSplit);
end

% ------------------------------------------------------------------------------
% Set project variables
% ------------------------------------------------------------------------------
parentdir = '~/Dropbox/Work/ResProjects/';
% parentdir = '/home/lindenmp/kg98_scratch/Linden/ResProjects/';
switch WhichProject
	case 'OCDPG'
		projdir = [parentdir,'rfMRI_denoise/OCDPG/'];
		sublist = [projdir,'OCDPGe.csv'];
		datadir = [projdir,'data/'];
        preprostr = '/func/prepro/';

		TR = 2.5;
	case 'UCLA'
		projdir = [parentdir,'rfMRI_denoise/UCLA/'];
		sublist = [projdir,'UCLA.csv'];
		datadir = [projdir,'data/'];
        preprostr = '/func/prepro/';

		TR = 2;
	case 'NYU_2'
		projdir = [parentdir,'rfMRI_denoise/NYU_2/'];
		sublist = [projdir,'NYU_2.csv'];
		datadir = [projdir,'data/'];
		% Baseline data directory string
		% Note, we use the baseline data to calculate motion
		preprostr = '/session_1/func_1/prepro/';
		% preprostr = '/session_1/func_2/prepro/';
		% preprostr = '/session_2/func_1/prepro/';
	
		TR = 2;
	case 'COBRE'
		projdir = [parentdir,'rfMRI_denoise/COBRE/'];
		sublist = [projdir,'COBRE.csv'];
		datadir = [projdir,'data/'];
        preprostr = '/func/prepro/';

		TR = 2;
	case 'Beijing_Zang'
		projdir = [parentdir,'rfMRI_denoise/Beijing_Zang/'];
		sublist = [projdir,'Beijing_Zang.csv'];
		datadir = [projdir,'data/'];
        preprostr = '/func/prepro/';

		TR = 2;
	case 'COBRE_HCTSA'
		projdir = [parentdir,'SCZ_HCTSA/COBRE/'];
		sublist = [projdir,'COBRE.csv'];
		datadir = [projdir,'data/'];
		preprostr = '/session_1/rest_1/prepro/';
	
		TR = 2;
	case 'UCLA_HCTSA'
		projdir = [parentdir,'SCZ_HCTSA/UCLA/'];
		sublist = [projdir,'UCLA.csv'];
		datadir = [projdir,'data/'];
		preprostr = '/func/prepro/';
	
		TR = 2;
	case 'NAMIC_HCTSA'
		projdir = [parentdir,'SCZ_HCTSA/NAMIC/'];
		sublist = [projdir,'NAMIC.csv'];
		datadir = [projdir,'data/'];
		preprostr = '/func/prepro/';
	
		TR = 3;
	case 'ML_SNS'
		projdir = '/projects/kg98/Michelle/PROJECTS/SNS/'
		sublist = [projdir,'NAMIC.csv'];
		datadir = projdir;
		preprostr = '/func/prepro/';

		TR = 2;
end

% ------------------------------------------------------------------------------
% Set parcellation
% Note, code is not setup to process multiple parcellations concurrently.
% ------------------------------------------------------------------------------
ROIDir = '~/Dropbox/Work/ROIs/';
% ROIDir = '/projects/kg98/Linden/ROIs/';
switch WhichParc
	case 'Gordon'
		Parc = 1;
		ROI_Coords = dlmread([ROIDir,'Gordon/Gordon_Centroids.txt']);
		fileName = [ROIDir,'Gordon/CommunityModified.txt'];
	case 'Power'
		Parc = 2;
		ROI_Coords = dlmread([ROIDir,'Power/Power2011_xyz_MNI.txt']);
		fileName = [ROIDir,'Power/Community.txt'];
end

fileID = fopen(fileName);
ROIStruct = textscan(fileID,'%s'); ROIStruct = ROIStruct{1};
% rearrange by community
[ROIStruct_com,ROI_idx] = sort(ROIStruct);

% Convert text labels to unique integer values
ROILabels = unique(ROIStruct);
ROIStructID = zeros(size(ROIStruct));
numROIComms = length(ROILabels);
for i = 1:numROIComms
	x = find(strcmp(ROIStruct, ROILabels{i}));
	ROIStructID(x) = i;
end


% ------------------------------------------------------------------------------
% Load ROI coordinates
% ------------------------------------------------------------------------------
% Calculate pairwise euclidean distance
ROIDist = pdist2(ROI_Coords,ROI_Coords,'euclidean');

% Flatten distance matrix
ROIDistVec = LP_FlatMat(ROIDist);

% Calculate number of ROIs
numROIs = size(ROIDist,1);

% Calculate number of edges
numConnections = numROIs * (numROIs - 1) / 2;

% ------------------------------------------------------------------------------
% 							Preprocessing pipelines
% ------------------------------------------------------------------------------
switch WhichProject
	case {'M3_COBRE','M3_UCLA','M3_NAMIC'}
		noiseOptions = {'sICA-AROMA+2P',...
						'sICA-AROMA+2P+GSR'};
		noiseOptionsNames = {'ICA-AROMA+2Phys',...
							'ICA-AROMA+2Phys+GSR'};
	case 'ML_SNS'
		noiseOptions = {'sICA-AROMA+2P'};
		noiseOptionsNames = {'ICA-AROMA+2Phys'};			
	otherwise
		noiseOptions = {'6P',...
						'6P+2P',...
						'6P+2P+GSR',...
						'24P',...
						'24P+8P',...
						'24P+8P+4GSR',...
						'24P+aCC',...
						'24P+aCC+4GSR',...
						'24P+aCC50',...
						'24P+aCC50+4GSR',...
						'12P+aCC',...
						'12P+aCC50',...
						'ICA-AROMA+2P',...
						'ICA-AROMA+2P+GSR',...
						'ICA-AROMA+8P',...
						'ICA-AROMA+8P+4GSR',...
						'24P+8P+4GSR+SpikeReg',...
						'24P+8P+4GSR+JP12Scrub',...
						'24P+4P+2GSR+JP14Scrub'};

		noiseOptionsNames = {'6HMP',...
							'6HMP+2Phys',...
							'6HMP+2Phys+GSR',...
							'24HMP',...
							'24HMP+8Phys',...
							'24HMP+8Phys+4GSR',...
							'24HMP+aCompCor',...
							'24HMP+aCompCor+4GSR',...
							'24HMP+aCompCor50',...
							'24HMP+aCompCor50+4GSR',...
							'12HMP+aCompCor',...
							'12HMP+aCompCor50',...
							'ICA-AROMA+2Phys',...
							'ICA-AROMA+2Phys+GSR',...
							'ICA-AROMA+8Phys',...
							'ICA-AROMA+8Phys+4GSR',...
							'24HMP+8Phys+4GSR+SpikeReg',...
							'24HMP+8Phys+4GSR+JP12Scrub',...
							'24HMP+4Phys+2GSR+JP14Scrub'};
end

numPrePro = length(noiseOptions);

% ------------------------------------------------------------------------------
% Subject list
% ------------------------------------------------------------------------------
fprintf(1, 'Loading metadata...\n');
metadata = readtable(sublist);
% convert participant IDs to strings
if ~iscellstr(metadata.ParticipantID)
	metadata.ParticipantID =  cellfun(@num2str, num2cell(metadata.ParticipantID), 'UniformOutput', false);
	if ismember('M3_COBRE',WhichProject,'rows')
		metadata.ParticipantID = cellfun(@(c)['00' c],metadata.ParticipantID,'UniformOutput',false);
	end
end

% Retain only group 1 (assumed to be HCs) and group 2 (assumed to be patients)
% I do this because the OCDPG dataset has some PGs in it that need to be removed
if excludeGroup
	% Only do this if there is more than a single group
	if numel(unique(metadata.Diagnosis)) > 1
		metadata(metadata.Diagnosis == 3,:) = [];
	end
end

numGroups = numel(unique(metadata.Diagnosis));

% ------------------------------------------------------------------------------
% Jenkinson's mean FD
% ------------------------------------------------------------------------------
fprintf(1, 'Loading Jenkinson''s mean FD metric\n');
[metadata.exclude,~,metadata.fdJenk,metadata.fdJenk_m] = GetExcludeForSample(datadir,metadata.ParticipantID,preprostr);
fprintf(1, 'done\n');

% compute number of volumes using the length of fdJenk
% note, this is assumed to be same for all subjects!
numVols = length(metadata.fdJenk{1});

% ------------------------------------------------------------------------------
% Perform initial exclusion based on gross movement
% ------------------------------------------------------------------------------
WhichExclude = 1;
fprintf(1, 'Excluded %u subjects based on gross movement\n', sum(metadata.exclude(:,WhichExclude)));
metadata(metadata.exclude(:,WhichExclude),:) = [];

% compute numsubs
numSubs = size(metadata,1);

% ------------------------------------------------------------------------------
% Plot Movement params
% ------------------------------------------------------------------------------
if numGroups > 1
	[h,p,~,stats] = ttest2(metadata.fdJenk_m(metadata.Diagnosis == 1),metadata.fdJenk_m(metadata.Diagnosis == 2),'Vartype','unequal');
	if p < 0.05
		fprintf(1, 'There is a significant group difference in mean FD. t-value = %s. p-value = %s\n',num2str(stats.tstat),num2str(p));
	end
end

BF_JitteredParallelScatter({metadata.fdJenk_m(metadata.Diagnosis == 1),metadata.fdJenk_m(metadata.Diagnosis == 2)});
ax = gca;
ax.XTick = [1,2];
ax.XTickLabel = {'Controls','Patients'};
xlabel('Group')
ylabel('Mean FD')

% ------------------------------------------------------------------------------
% Variables
% ------------------------------------------------------------------------------
	allData = struct('noiseOptions',noiseOptions,...
					'noiseOptionsNames',noiseOptionsNames,...
					'cfg',struct,...
					'FC',[],...
					'FCVec',[],...
					'VarCovar',[],...
					'Var',[],...
					'GCOR',[],...
					'NaNFilter',[],...
					'QCFCVec',[],...
					'QCFC_PropSig_corr',[],...
					'QCFC_PropSig_unc',[],...
					'QCFC_DistDep',[],...
					'QCFC_DistDep_Pval',[],...
					'MeanEdgeWeight',[],...
					'tDOF',[],...
					'tDOF_mean',[],...
					'tDOF_std',[],...
					'tDOF_gmean',[],...
					'tDOF_gstd',[]);

	% Add censoring variables if runScrub or runSR are true
	if runScrub | runSR
		% censored FC
		allData().FC_censored = [];

		% delta r
		allData().FC_delta_mean = [];

		% pre-censoring
		allData().preCensor_QCFCVec = [];
		allData().preCensor_QCFC_DistDep = [];
		allData().preCensor_QCFC_DistDep_Pval = [];
		allData().preCensor_QCFC_PropSig_corr = [];
		allData().preCensor_QCFC_PropSig_unc = [];

		% post-censoring
		allData().postCensor_QCFCVec = [];
		allData().postCensor_QCFC_DistDep = [];
		allData().postCensor_QCFC_DistDep_Pval = [];
		allData().postCensor_QCFC_PropSig_corr = [];
		allData().postCensor_QCFC_PropSig_unc = [];
	end

	if ismember('OCDPG',WhichProject,'rows') | ismember('UCLA',WhichProject,'rows')
		allData().NBS_statMat = [];
		allData().NBS_statVec = [];

		allData().NBS_sigMat = [];
		allData().NBS_sigVec = [];

		allData().NaNFilterMat = [];

		allData().NBS_PropSig = [];
		
		allData().u05_PropSig = [];
		allData().u001_PropSig = [];
		allData().FDR_PropSig = [];
	end

	% Add test-retest variables if WhichProject == 'NYU_2'
	if ismember('NYU_2',WhichProject,'rows')
		allData().cfg_wrt = struct;
		allData().cfg_brt = struct;
		allData().FCw = [];
		allData().FCwVec = [];
		allData().FCb = [];
		allData().FCbVec = [];
		allData().ICCw = [];
		allData().ICCw_mean = [];
		allData().ICCw_std = [];
		allData().ICCb = [];
		allData().ICCb_mean = [];
		allData().ICCb_std = [];
		allData().ICC_tstat = [];
		allData().ICC_Pval = [];
		allData().ICC_Pval_corr = [];
	end

% ------------------------------------------------------------------------------
% Power's scrubbing
% ------------------------------------------------------------------------------
if runScrub
	fprintf(1, 'Loading scrubbing mask\n');
	[metadata.ScrubMask,metadata.mov,metadata.fdPower,metadata.dvars] = GetScrubbingForSample(datadir,metadata.ParticipantID,preprostr);
	fprintf(1, 'done\n');
end

% ------------------------------------------------------------------------------
% Create exclusion vector based on <4 minutes of post-scrub BOLD data
% Note, subjects not excluded yet. that is done during the preprocessing loop below
% ------------------------------------------------------------------------------
if runScrub | runSR
	% exclusion vector
	excludeCensor = zeros(numSubs,1);
	% threshold for exclusion in minutes
	thresh = 4;

	% Loop over subjects and check if censoring would reduce data to <4 minutes
	for i = 1:numSubs
		% Find number of vols left after censoring
		if runScrub
			% number of volumes - number of scrubbed volumes
			numCVols = numVols - sum(metadata.ScrubMask{i});			
		elseif runSR
			% Get spike regressors for subject i
			spikereg = GetSpikeRegressors(metadata.fdJenk{i},0.25);
			% number of volumes - number of spike regressors (columns)
			numCVols = numVols - size(spikereg,2);
		end	

		% Compute length, in minutes, of time series data left after censoring
		NTime = (numCVols * TR)/60;
		% if less than threshold, mark for exclusion
		if NTime < thresh
			excludeCensor(i) = 1;
		end
	end
	% convert to logical
	excludeCensor = logical(excludeCensor);

	fprintf(1, '%u subjects marked for exclusion based on 4 minute censoring criteria \n', sum(excludeCensor));
end

% ------------------------------------------------------------------------------
% 						Loop over preprocessing pipelines
% ------------------------------------------------------------------------------
for i = 1:numPrePro
    removeNoise = allData(i).noiseOptions;
    removeNoiseName = allData(i).noiseOptionsNames;
	fprintf(1, '\nProcessing data: %s\n',removeNoise);

	% ------------------------------------------------------------------------------
	% Get time series and functional connectivity data
	% ------------------------------------------------------------------------------
	cfgFile = 'cfg.mat';
	[allData(i).cfg,allData(i).FC,allData(i).FCVec,allData(i).VarCovar,allData(i).Var,allData(i).GCOR] = GetFCForSample(datadir,metadata.ParticipantID,preprostr,removeNoise,cfgFile,Parc,numROIs,numConnections);

	% ------------------------------------------------------------------------------
	% Compute QC-FC
	% Note, make sure use z transformed r values
	% ------------------------------------------------------------------------------
	fprintf(1, 'Computing QC-FC: %s\n',removeNoise);
	[allData(i).QCFCVec,allData(i).NaNFilter,allData(i).QCFC_PropSig_corr,allData(i).QCFC_PropSig_unc,~,allData(i).QCFC_DistDep,allData(i).QCFC_DistDep_Pval] = RunQCFC(metadata.fdJenk_m,allData(i).FC,ROIDistVec);

	% ------------------------------------------------------------------------------
	% Compute mean edge weight
	% ------------------------------------------------------------------------------
	fprintf(1, 'Computing mean edge weights: %s\n',removeNoise);
	FCVec = [];
	for j = 1:numSubs
		% flatten FC for subject j
		vec = LP_FlatMat(allData(i).FC(:,:,j));
		% filter NaNs from QCFC analyses
		vec = vec(allData(i).NaNFilter);
		% store
		FCVec(:,j) = vec;
	end

	% average across subject
	allData(i).MeanEdgeWeight = mean(FCVec,2);

	% ------------------------------------------------------------------------------
	% Censoring
	% ------------------------------------------------------------------------------
	if runScrub | runSR
		fprintf(1, 'Performing censoring analysis: %s\n',removeNoise);

	    if runSR
			for j = 1:numSubs
				% Load in spike regressed time series data
			    clear temp
			    temp = load([datadir,metadata.ParticipantID{j},preprostr,removeNoise,'+SpikeReg/cfg.mat']);
			    
			    allData(i).cfg(j).roiTS_censored{Parc} = temp.cfg.roiTS{Parc};
			    allData(i).cfg(j).noiseTS_spikereg = temp.cfg.noiseTS;
			end
		end

		if runScrub
			% scrub time series
			for j = 1:numSubs
				allData(i).cfg(j).roiTS_censored{Parc} = allData(i).cfg(j).roiTS{Parc}(~metadata.ScrubMask{j},:);
			end
		end

		% compute correlations
		allData(i).FC_censored = zeros(numROIs,numROIs,numSubs); 
		for j = 1:numSubs
			allData(i).FC_censored(:,:,j) = corr(allData(i).cfg(j).roiTS_censored{Parc});
			% Perform fisher z transform
			allData(i).FC_censored(:,:,j) = atanh(allData(i).FC_censored(:,:,j));
		end

		% ------------------------------------------------------------------------------
		% Compute delta correlations
		% ------------------------------------------------------------------------------
		% find difference between correlations before and after scrubbing
		FC_delta = allData(i).FC_censored(:,:,~excludeCensor) - allData(i).FC(:,:,~excludeCensor);

		% Take average delta
		allData(i).FC_delta_mean = nanmean(FC_delta,3);

		% ------------------------------------------------------------------------------
		% Compute pre-censoring QC-FC
		% ------------------------------------------------------------------------------
		[allData(i).preCensor_QCFCVec,~,allData(i).preCensor_QCFC_PropSig_corr,allData(i).preCensor_QCFC_PropSig_unc,~,allData(i).preCensor_QCFC_DistDep,allData(i).preCensor_QCFC_DistDep_Pval] = RunQCFC(metadata.fdJenk_m(~excludeCensor),allData(i).FC(:,:,~excludeCensor),ROIDistVec);

		% ------------------------------------------------------------------------------
		% Compute post-censoring QC-FC
		% ------------------------------------------------------------------------------
		[allData(i).postCensor_QCFCVec,~,allData(i).postCensor_QCFC_PropSig_corr,allData(i).postCensor_QCFC_PropSig_unc,~,allData(i).postCensor_QCFC_DistDep,allData(i).postCensor_QCFC_DistDep_Pval] = RunQCFC(metadata.fdJenk_m(~excludeCensor),allData(i).FC_censored(:,:,~excludeCensor),ROIDistVec);
	end

	% ------------------------------------------------------------------------------
	% Get tDOF
	% ------------------------------------------------------------------------------
	fprintf(1, 'Computing tDOF: %s\n',removeNoise);
	allData(i).tDOF = zeros(numSubs,1);
	for j = 1:numSubs

		% get tDOF
		% First, find size of second dimension of noiseTS
		if runSR
			allData(i).tDOF(j) = size(allData(i).cfg(j).noiseTS_spikereg,2);
		else
			allData(i).tDOF(j) = size(allData(i).cfg(j).noiseTS,2);
		end

		if runScrub
			allData(i).tDOF(j) = allData(i).tDOF(j) + sum(metadata.ScrubMask{j});
		end

		% Then, if ICA-AROMA pipeline, find number of ICs and add to tDOF
		if ~isempty(strfind(removeNoise,'ICA-AROMA'))
			if runSR
				x = dlmread([datadir,metadata.ParticipantID{j},preprostr,removeNoise,'+SpikeReg/classified_motion_ICs.txt']);
			else
				x = dlmread([datadir,metadata.ParticipantID{j},preprostr,removeNoise,'/classified_motion_ICs.txt']);
			end
			allData(i).tDOF(j) = allData(i).tDOF(j) + length(x);
		end
	end

	% ------------------------------------------------------------------------------
	% Calculate mean temporal degrees of freedom lost
	% ------------------------------------------------------------------------------
	% tDOF will be the same for most pipelines
	% but some have variable regressor amounts
	% so we take mean over subjects
	if runScrub | runSR
		allData(i).tDOF_mean = mean(allData(i).tDOF(~excludeCensor));
		allData(i).tDOF_std = std(allData(i).tDOF(~excludeCensor));
	else
		allData(i).tDOF_mean = mean(allData(i).tDOF);
		allData(i).tDOF_std = std(allData(i).tDOF);
	end

	% also get mean by diagnostic groups (metadata.Diagnosis)
	if ismember('OCDPG',WhichProject,'rows') | ismember('UCLA',WhichProject,'rows')
		allData(i).tDOF_gmean(1) = mean(allData(i).tDOF(metadata.Diagnosis == 1));
		allData(i).tDOF_gmean(2) = mean(allData(i).tDOF(metadata.Diagnosis == 2));
		
		allData(i).tDOF_gstd(1) = std(allData(i).tDOF(metadata.Diagnosis == 1));
		allData(i).tDOF_gstd(2) = std(allData(i).tDOF(metadata.Diagnosis == 2));
	end

	% ------------------------------------------------------------------------------
	% Perform t-test on tDOF-loss
	% ------------------------------------------------------------------------------
	if ismember('OCDPG',WhichProject,'rows') | ismember('UCLA',WhichProject,'rows')
		x = allData(i).tDOF(metadata.Diagnosis == 1);
		y = allData(i).tDOF(metadata.Diagnosis == 2);

		[h,p,~,stats] = ttest2(x,y,'Vartype','unequal');
		if p < 0.05
			fprintf(1, 'Significant group difference in tDOF-loss. t-value = %s. p-value = %s\n',num2str(stats.tstat),num2str(p));
		else
			fprintf(1, 'NO significant group difference in tDOF-loss. t-value = %s. p-value = %s\n',num2str(stats.tstat),num2str(p));
		end
		fprintf(1, '\tMean tDOF-loss, group 1: %s\n', num2str(round(allData(i).tDOF_gmean(1),2)));
		fprintf(1, '\tMean tDOF-loss, group 2: %s\n', num2str(round(allData(i).tDOF_gmean(2),2)));
	end

	% ------------------------------------------------------------------------------
	% Get NBS contrasts
	% ------------------------------------------------------------------------------
	if ismember('OCDPG',WhichProject,'rows') | ismember('UCLA',WhichProject,'rows')
		fprintf(1, 'Getting NBS contrasts: %s\n',removeNoise);
	    % ------------------------------------------------------------------------------
	    % Load NBS data
	    % ------------------------------------------------------------------------------
	    Tval = 1; % only 1 works for now.

	    file = dir([nbsdir,WhichSplit,'_',WhichParc,'_',removeNoise,'_NBS_t',num2str(Tval),'*.mat']);
	    load([nbsdir,file(1).name])

	    % ------------------------------------------------------------------------------
	    % Loop over NBS contrasts
	    % ------------------------------------------------------------------------------
	    numContrasts = 2;
	    for j = 1:numContrasts
		    % Get nbs data for jth contrast
		    nbs = nbsOut{j};
		    
		    % ------------------------------------------------------------------------------
		    % Initialise
		    % ------------------------------------------------------------------------------
			% Matrix of test stats
			allData(i).NBS_statMat{j} = zeros(numROIs,numROIs);
			% vector of above
			allData(i).NBS_statVec{j} = zeros(numConnections,1);

			% Binary of sig edges in statMat
			allData(i).NBS_sigMat{j} = zeros(numROIs,numROIs);
			% vector of above
			allData(i).NBS_sigVec{j} = zeros(numConnections,1);

			% ------------------------------------------------------------------------------
			% Get stat Mat
			% ------------------------------------------------------------------------------
			fprintf(1, '\tFound %u networks for contrast %u \n', length(nbs.NBS.con_mat),j);

		    % matrix of test statistics
			allData(i).NBS_statMat{j} = nbs.NBS.test_stat;

			% flatten upper triangle
			% vector of unthrehsolded test statistics for distance plots
			allData(i).NBS_statVec{j} = LP_FlatMat(allData(i).NBS_statMat{j}); 

			% ------------------------------------------------------------------------------
			% Convert NanFilter to matrix mask
			% ------------------------------------------------------------------------------
			x = LP_SquareVec(allData(i).NaNFilter,numROIs);
			allData(i).NaNFilterMat = [x+x'];
			y = sum(allData(i).NaNFilterMat); y(y > 0) = 1;
			allData(i).NaNFilterMat(eye(numROIs) == 1) = y;
			allData(i).NaNFilterMat = logical(allData(i).NaNFilterMat);

			% ------------------------------------------------------------------------------
			% Threshold statMat for matrix plot
			% ------------------------------------------------------------------------------
			if ~isempty(nbs.NBS.con_mat)
				% binary matrix denoting significant edges
				allData(i).NBS_sigMat{j} = nbs.NBS.con_mat{1};
				% get bottom triangle for accurate reordering of rows and columns (NBS only outputs upper triangle)
				allData(i).NBS_sigMat{j} = [allData(i).NBS_sigMat{j} + allData(i).NBS_sigMat{j}'];
				% flatten upper triangle
				allData(i).NBS_sigVec{j} = LP_FlatMat(allData(i).NBS_sigMat{j});
			end

			% ------------------------------------------------------------------------------
			% Calculate proportion of significant connections adjusted by NaNFilter
			% ------------------------------------------------------------------------------
		    allData(i).NBS_PropSig(j) = sum(allData(i).NBS_sigVec{j}(allData(i).NaNFilter)) / sum(allData(i).NaNFilter) * 100;
		    % allData(i).NBS_PropSig(j) = sum(allData(i).NBS_sigVec{j}) / numConnections * 100;

			% ------------------------------------------------------------------------------
			% Calculate PropSig at different thresholds for comparison
			% ------------------------------------------------------------------------------
			% Get design matrix w/ Cov from NBS
		    X = nbs.GLM.X;

		    % number of subjects in NBS
		    numSubsNBS = size(X,1);

		    % Get NBS t stat
			t_stat = LP_FlatMat(allData(i).NBS_statMat{j});
			% Calculate p-value
			pval = 1 - tcdf(t_stat,numSubsNBS - rank(X));

			% 1) p < 0.05 unc
			pval_thr = pval < 0.05;
		    allData(i).u05_PropSig(j) = sum(pval_thr(allData(i).NaNFilter)) / sum(allData(i).NaNFilter) * 100;

			% 2) p < 0.001 unc
			pval_thr = pval < 0.001;
		    allData(i).u001_PropSig(j) = sum(pval_thr(allData(i).NaNFilter)) / sum(allData(i).NaNFilter) * 100;

			% 3) p < 0.05 FDR
			pval_cor = mafdr(pval,'BHFDR','true');
			pval_thr = pval_cor < 0.05;
		    allData(i).FDR_PropSig(j) = sum(pval_thr(allData(i).NaNFilter)) / sum(allData(i).NaNFilter) * 100;
		end
		% convert PropSig to full double. (for some reason it comes sparse)
		allData(i).NBS_PropSig = full(allData(i).NBS_PropSig);

		% clean up
		clear nbsOut nbs
	end

	% ------------------------------------------------------------------------------
	% Compute ICC
	% ------------------------------------------------------------------------------
	if ismember('NYU_2',WhichProject,'rows')
		% ------------------------------------------------------------------------------
		% Within session retest data
		% Get time series and functional connectivity data
		% ------------------------------------------------------------------------------
		[allData(i).cfg_wrt,allData(i).FCw,allData(i).FCwVec,~,~] = GetFCForSample(datadir,ParticipantIDs,'/session_1/rest_2/prepro/',removeNoise,cfgFile,Parc,numROIs,numConnections);

		% ------------------------------------------------------------------------------
		% Between session retest data
		% Get time series and functional connectivity data
		% ------------------------------------------------------------------------------
		[allData(i).cfg_brt,allData(i).FCb,allData(i).FCbVec,~,~] = GetFCForSample(datadir,ParticipantIDs,'/session_2/rest_1/prepro/',removeNoise,cfgFile,Parc,numROIs,numConnections);

		% ------------------------------------------------------------------------------
		% Calculate ICC
		% ------------------------------------------------------------------------------
		fprintf(1, 'Performing ICC for: %s',removeNoise);

		% initialise
		allData(i).ICCw = zeros(1,numConnections);
		allData(i).ICCb = zeros(1,numConnections);

		% For each functional connection
		for j = 1:numConnections
			% 1) within session ICC
			x = [allData(i).FCVec(:,j), allData(i).FCwVec(:,j)];
			allData(i).ICCw(j) = GetICC(x);

			% 2) between session ICC
			y = [allData(i).FCVec(:,j), allData(i).FCbVec(:,j)];
			allData(i).ICCb(j) = GetICC(y);
		end

		% filter nans
		allData(i).ICCw = allData(i).ICCw(allData(i).NaNFilter);
		allData(i).ICCb = allData(i).ICCb(allData(i).NaNFilter);

		% mean
		allData(i).ICCw_mean = nanmean(allData(i).ICCw);
		allData(i).ICCw_std = nanstd(allData(i).ICCw);

		% std
		allData(i).ICCb_mean = nanmean(allData(i).ICCb);
		allData(i).ICCb_std = nanstd(allData(i).ICCb);

		fprintf(1, '...done \n');

		% ------------------------------------------------------------------------------
		% Compute within session and between session ICC differences
		% ------------------------------------------------------------------------------
		WhichStat = 't';
		switch WhichStat
			case 't'
				[~,p,~,stats] = ttest(allData(i).ICCw,allData(i).ICCb);

				allData(i).ICC_tstat = stats.tstat;
				allData(i).ICC_Pval = p;
		end
	end
end

% ------------------------------------------------------------------------------
% Check whether NaNs filtered out of each pipeline are the same.
% They should be since if there are NaNs in QCFC, the same subject(s) will cause them each time.
% But, still worth checking explicitly
% ------------------------------------------------------------------------------
temp = [allData(:).NaNFilter]';
allRowsEqual = size(unique(temp,'rows'),1) == 1;
if allRowsEqual == 1
	% If all nan filters are the same across pipelines
	% Safe to use any of the nan filters to filter ROI dist vec permenantly
    % ROIDistVec = ROIDistVec(allData(end).NaNFilter);
elseif allRowsEqual ~= 1
    warning('NaN filters are not the same across pipelines!');
    % ROIDistVec = ROIDistVec(allData(end).NaNFilter);
end

% ------------------------------------------------------------------------------
% Correct ICC t-stats for multiple comparisons
% ------------------------------------------------------------------------------
if ismember('NYU_2',WhichProject,'rows')
	% correct p values using FDR
	x = mafdr([allData(:).ICC_Pval],'BHFDR','true');
	for i = 1:numPrePro
		allData(i).ICC_Pval_corr = x(i);
	end
end

% ------------------------------------------------------------------------------
% Figures
% ------------------------------------------------------------------------------
FSize = 10;

if runPlot
	clear extraParams
	% ------------------------------------------------------------------------------
	% Chart colors and line styles
	% ------------------------------------------------------------------------------
	tempColors = num2cell([255,105,97;97,168,255;178,223,138;117,112,179]./255,2);  
	for i = 1:numPrePro
		strs = strsplit(allData(i).noiseOptions,'+');
        % colors
        if any(strmatch('6P',strs,'exact')) == 1; theColors{i} = tempColors{1}; end
        if any(strmatch('24P',strs,'exact')) == 1; theColors{i} = tempColors{2}; end
        if any(strmatch('aCC',strs,'exact')) == 1 || any(strmatch('aCC50',strs,'exact')) == 1; theColors{i} = tempColors{3}; end
        if any(strmatch('sICA-AROMA',strs,'exact')) == 1; theColors{i} = tempColors{4}; end
        
        % line style
        if any(strmatch('GSR',strs,'exact')) == 1 || any(strmatch('4GSR',strs,'exact')) == 1
        	theLines{i} = ':';
        	% theLines{i} = '--';
        else
        	theLines{i} = '-';
        end
	end

	% ------------------------------------------------------------------------------
	% QC-FC significant proportion
	% Corresponds to Figure 1 and Figure 4 in case of censoring
	% ------------------------------------------------------------------------------
	if ~runScrub & ~runSR
		Fig_QCFC_Dist = figure('color','w', 'units', 'centimeters', 'pos', [0 0 16 9], 'name',['Fig_QCFC_Dist']); box('on'); movegui(Fig_QCFC_Dist,'center');
		sp = subplot(1,3,1);
		pos = get(sp,'Position');
	    set(gca,'Position',[pos(1)*2.25, pos(2)*1.2, pos(3)*1.5, pos(4)*1]); % [left bottom width height]

		% Create data
		data = {[allData(:).QCFC_PropSig_corr]'};
		data_std = cell(1,length(data)); [data_std{:}] = deal(zeros(size(data{1})));

		% Create table
		T = table(data{1},'RowNames',{allData(:).noiseOptionsNames}','VariableNames',{'QCFC_PropSig'})

		% Create bar chart
		clear extraParams
		extraParams.xTickLabels = {allData(:).noiseOptionsNames};
		extraParams.xLabel = ''; % 'Pipeline'
		extraParams.yLabel = 'QC-FC (%)';
		extraParams.plotWidth = 10.75;
		extraParams.plotHeight = 10;
	    extraParams.theColors = theColors;
	    extraParams.theLines = theLines;
		extraParams.yLimits = [0 100];

		TheBarChart(data,data_std,false,extraParams)

		% ------------------------------------------------------------------------------
		% QCFC distributions
		% ------------------------------------------------------------------------------
		sp = subplot(1,3,3);
		pos = get(sp,'Position');
	    set(gca,'Position',[pos(1)*0.925, pos(2)*1.2, pos(3)*1.55, pos(4)*1]); % [left bottom width height]
		clear extraParams
		% extraParams.theLabels = {allData(:).noiseOptionsNames};
		extraParams.customSpot = '';
		extraParams.add0Line = true;
	    extraParams.theColors = theColors;
		BF_JitteredParallelScatter({allData(:).QCFCVec},1,1,0,extraParams)
		ax = gca;

		% Set axis stuff
		ax.FontSize = FSize;
		% XTicks = [1:size(data{1},1)];
		XTicks = [];
		ax.XTick = XTicks;
		ax.XTickLabel = [];
		ax.XLim = ([0 numPrePro+1]);
		ax.YLim = ([-0.75 0.75]);
		ylabel('QC-FC')

		view(90,90)
	elseif runScrub | runSR
		% Create data
		data = {[allData(:).QCFC_PropSig_corr;allData(:).preCensor_QCFC_PropSig_corr;allData(:).postCensor_QCFC_PropSig_corr]'};
		data_std = cell(1,length(data)); [data_std{:}] = deal(zeros(size(data{1})));

		% Create table
		T = table(data{1},'RowNames',{allData(:).noiseOptionsNames}','VariableNames',{'QCFC_PropSig'})

		% Create bar chart
		clear extraParams
		% if runScrub
		% 	extraParams.Legend = {'Before exclusion','Pre-Scrubbing','Post-Scrubbing'};
		% elseif runSR
		% 	extraParams.Legend = {'Before exclusion','Pre-SpikeReg','Post-SpikeReg'};
		% end
		extraParams.xTickLabels = {allData(:).noiseOptionsNames};
		extraParams.xLabel = 'Pipeline';
		extraParams.yLabel = 'QC-FC (%)';
		extraParams.plotWidth = 10.5;
		extraParams.plotHeight = 9;
		extraParams.yLimits = [0 25];

		TheBarChart(data,data_std,true,extraParams)
	end

	% ------------------------------------------------------------------------------
	% QC-FC uncorrected
	% ------------------------------------------------------------------------------
	if ~runScrub & ~runSR
		% Create data
		data = {[allData(:).QCFC_PropSig_unc]'};
		data_std = cell(1,length(data)); [data_std{:}] = deal(zeros(size(data{1})));

		% Create table
		T = table(data{1},'RowNames',{allData(:).noiseOptionsNames}','VariableNames',{'QCFC_PropSig_unc'})

		% Create bar chart
		clear extraParams
		extraParams.xTickLabels = {allData(:).noiseOptionsNames};
		extraParams.xLabel = ''; % 'Pipeline'
		extraParams.yLabel = 'QC-FC (%)';
		extraParams.plotWidth = 10.5;
		extraParams.plotHeight = 9;
	    extraParams.theColors = theColors;
	    extraParams.theLines = theLines;
		extraParams.yLimits = [0 100];

		TheBarChart(data,data_std,true,extraParams)
	elseif runScrub | runSR
		% Create data
		data = {[allData(:).QCFC_PropSig_unc;allData(:).preCensor_QCFC_PropSig_unc;allData(:).postCensor_QCFC_PropSig_unc]'};
		data_std = cell(1,length(data)); [data_std{:}] = deal(zeros(size(data{1})));

		% Create table
		T = table(data{1},'RowNames',{allData(:).noiseOptionsNames}','VariableNames',{'QCFC_PropSig_unc'})

		% Create bar chart
		clear extraParams
		extraParams.xTickLabels = {allData(:).noiseOptionsNames};
		extraParams.xLabel = 'Pipeline';
		extraParams.yLabel = 'QC-FC (%)';
		extraParams.plotWidth = 10.5;
		extraParams.plotHeight = 9;
		extraParams.yLimits = [0 50];

		TheBarChart(data,data_std,true,extraParams)
	end

	% ------------------------------------------------------------------------------
	% QC-FC distance dependence
	% Corresponds to Figure 2 and Figure 4 in case of censoring
	% ------------------------------------------------------------------------------
	% Create data
	if ~runScrub & ~runSR
		data = {[allData(:).QCFC_DistDep]'};
	elseif runScrub | runSR
		data = {[allData(:).QCFC_DistDep;allData(:).preCensor_QCFC_DistDep;allData(:).postCensor_QCFC_DistDep]'};
	end
	data_std = cell(1,length(data)); [data_std{:}] = deal(zeros(size(data{1})));
	
	% Create table
	T = table(data{1},'RowNames',{allData(:).noiseOptionsNames}','VariableNames',{'QCFC_DistDep'})

	% Create bar chart
	extraParams.xTickLabels = {allData(:).noiseOptionsNames};
	extraParams.xLabel = 'Pipeline';
	extraParams.yLabel = 'QC-FC distance dependence';
	extraParams.plotWidth = 10.5;
	extraParams.plotHeight = 9;
    extraParams.theColors = theColors;
    extraParams.theLines = theLines;
	extraParams.yLimits = [-0.5 0.5];

	TheBarChart(data,data_std,true,extraParams)

	% ------------------------------------------------------------------------------
	% tDOF
	% Corresponds to Figure 5
	% ------------------------------------------------------------------------------
	% Create data
	data = {[allData(:).tDOF_mean]'};
	data_std = {[allData(:).tDOF_std]'};

	% Create table
	T = table(data{1},'RowNames',{allData(:).noiseOptionsNames}','VariableNames',{'tDOF_mean'})

	% Create bar chart
	extraParams.xTickLabels = {allData(:).noiseOptionsNames};
	extraParams.xLabel = 'Pipeline';
	extraParams.yLabel = 'tDOF-loss';
	extraParams.plotWidth = 10.5;
	extraParams.plotHeight = 9;
    extraParams.theColors = theColors;
    extraParams.theLines = theLines;
	extraParams.yLimits = [0 165];

	TheBarChart(data,data_std,true,extraParams)

	% ------------------------------------------------------------------------------
	% Size of significant edge component
	% Corresponds to Figure 6
	% ------------------------------------------------------------------------------
	if ismember('OCDPG',WhichProject,'rows') | ismember('UCLA',WhichProject,'rows')
		% Create table
	    if ismember('Diagnostic',WhichSplit,'rows')
			PropSig = vertcat(allData(:).NBS_PropSig);
			% PropSig = vertcat(allData(:).FDR_PropSig);
		elseif ismember('Motion',WhichSplit,'rows')
			PropSig = vertcat(allData(:).u05_PropSig);
			% PropSig = vertcat(allData(:).u001_PropSig);
			% PropSig = vertcat(allData(:).FDR_PropSig);
		end

		data = {PropSig(:,1),PropSig(:,2)};
		data_std = cell(1,length(data)); [data_std{:}] = deal(zeros(size(data{1})));
		T = table([data{:}],'RowNames',{allData(:).noiseOptionsNames}','VariableNames',{'PropSig'})

		% Create bar chart
		extraParams.xTickLabels = {allData(:).noiseOptionsNames};
		extraParams.xLabel = 'Pipeline';
		extraParams.yLabel = 'Significant connections (%)';
		extraParams.plotWidth = 10.5;
		extraParams.plotHeight = 9;
	    extraParams.theColors = theColors;
	    extraParams.theLines = theLines;
	    extraParams.makeABS = true;
	    if ismember('Diagnostic',WhichSplit,'rows')
			if ismember('OCDPG',WhichProject,'rows')
				extraParams.yLimits = [-20 20];
			elseif ismember('UCLA',WhichProject,'rows')
				extraParams.yLimits = [-45 45];
			end
		elseif ismember('Motion',WhichSplit,'rows')
			if ismember('OCDPG',WhichProject,'rows')
				extraParams.yLimits = [-20 20];
			elseif ismember('UCLA',WhichProject,'rows')
				extraParams.yLimits = [-100 40];
			end
	    end

		TheBarChart(data,data_std,true,extraParams)
	end

	% ------------------------------------------------------------------------------
	% ICC
	% Corresponds to Figure 7
	% ------------------------------------------------------------------------------
	if ismember('NYU_2',WhichProject,'rows')
		% ------------------------------------------------------------------------------
		% Sort tDOF
		% ------------------------------------------------------------------------------
		[~,tDOF_idx] = sort([allData(:).tDOF_mean],'ascend');

		% ------------------------------------------------------------------------------
		% WRT ICC
		% ------------------------------------------------------------------------------
		Fig_ICCw = figure('color','w', 'units', 'centimeters', 'pos', [0 0 10.5 9], 'name',['Fig_ICCw']); box('on'); movegui(Fig_ICCw,'center');
		data = {allData(:).ICCw}';
		RowNames = {allData(:).noiseOptionsNames}';

		clear extraParams
		extraParams.customSpot = '';
		extraParams.add0Line = true;
		extraParams.theColors = theColors(tDOF_idx);
		BF_JitteredParallelScatter(data(tDOF_idx),1,1,0,extraParams)
		ax = gca;

		% Set axis stuff
		ax.FontSize = FSize;
		ax.XTick = 1:numPrePro;
		ax.XTickLabel = RowNames(tDOF_idx);
		ax.XLim = ([0 numPrePro+1]);
		ax.YLim = ([-1 1]);
		xlabel('Pipeline')
		ylabel('Within Session ICC')
		view(90,90)

		% ------------------------------------------------------------------------------
		% BRT ICC
		% ------------------------------------------------------------------------------
		Fig_ICCb = figure('color','w', 'units', 'centimeters', 'pos', [0 0 10.5 9], 'name',['Fig_ICCb']); box('on'); movegui(Fig_ICCb,'center');
		data = {allData(:).ICCb}';
		RowNames = {allData(:).noiseOptionsNames}';

		clear extraParams
		extraParams.customSpot = '';
		extraParams.add0Line = true;
		extraParams.theColors = theColors(tDOF_idx);
		BF_JitteredParallelScatter(data(tDOF_idx),1,1,0,extraParams)
		ax = gca;

		% Set axis stuff
		ax.FontSize = FSize;
		ax.XTick = 1:numPrePro;
		ax.XTickLabel = RowNames(tDOF_idx);
		ax.XLim = ([0 numPrePro+1]);
		ax.YLim = ([-1 1]);
		xlabel('Pipeline')
		ylabel('Between Session ICC')
		view(90,90)

		% ------------------------------------------------------------------------------
		% ICC tstat
		% ------------------------------------------------------------------------------
		% Create data
		data = {[allData(:).ICC_tstat]'};
		data_std = cell(1,length(data)); [data_std{:}] = deal(zeros(size(data{1})));
		
		% reorder by tDOF-loss
		data{1} = data{1}(tDOF_idx);
		data_std{1} = data_std{1}(tDOF_idx);

		% Create table
		RowNames = {allData(:).noiseOptionsNames}';
		T = table(data{1},'RowNames',RowNames(tDOF_idx),'VariableNames',{'ICC_tstat'})

		% Create bar chart
		extraParams.plotWidth = 10.5;
		extraParams.plotHeight = 9;
		extraParams.xTickLabels = RowNames(tDOF_idx);
		extraParams.xLabel = 'Pipeline';
		extraParams.yLabel = 'ICC t-stat';
		extraParams.yLimits = [0 230];
	    extraParams.theColors = theColors(tDOF_idx);
    	extraParams.theLines = theLines;

		TheBarChart(data,data_std,true,extraParams)
	end
end

% ------------------------------------------------------------------------------
% Big Figures
% ------------------------------------------------------------------------------
if runBigPlots
	% ------------------------------------------------------------------------------
	% Figures: QCFC
	% ------------------------------------------------------------------------------
	% Initialise figures
	Fig_QCFC_DistDep = figure('color','w', 'units', 'centimeters', 'pos', [0 0 25 27], 'name',['Fig_QCFC_DistDep']); box('on'); movegui(Fig_QCFC_DistDep,'center');
	Fig_QCFC_MeanEdgeWeight = figure('color','w', 'units', 'centimeters', 'pos', [0 0 25 27], 'name',['Fig_QCFC_MeanEdgeWeight']); box('on'); movegui(Fig_QCFC_MeanEdgeWeight,'center');
	Fig_QCFC_tSNR = figure('color','w', 'units', 'centimeters', 'pos', [0 0 25 27], 'name',['Fig_QCFC_tSNR']); box('on'); movegui(Fig_QCFC_tSNR,'center');

	for i = 1:numPrePro
	    removeNoise = allData(i).noiseOptions;
	    removeNoiseName = allData(i).noiseOptionsNames;

		% ------------------------------------------------------------------------------
		% Plot: distance dependence
		% ------------------------------------------------------------------------------
		figure(Fig_QCFC_DistDep)
		subplot(4,ceil(numPrePro/4),i)
		set(gca,'FontSize',FSize)

		% Bin QCFC data by distance and generate means and stds for each
		numThresholds = 11;
		BF_PlotQuantiles(ROIDistVec(allData(i).NaNFilter),allData(i).QCFCVec,numThresholds)
		hold on
		plot([0:200],zeros(1,201),'--','Color','k')

		if ismember('OCDPG',WhichProject,'rows')
			ylim([-0.1 0.2])
		elseif ismember('UCLA',WhichProject,'rows')
			ylim([-0.1 0.3])
		elseif ismember('NYU_2',WhichProject,'rows')
			ylim([-0.1 0.3])
		end
		
		xlabel('Distance (mm)')
		ylabel('QC-FC')

		title(removeNoiseName,'Interpreter', 'none','FontSize',FSize,'FontWeight','normal')

		% ------------------------------------------------------------------------------
		% Plot: mean edge weight
		% ------------------------------------------------------------------------------
		figure(Fig_QCFC_MeanEdgeWeight)
		subplot(4,ceil(numPrePro/4),i)
		set(gca,'FontSize',FSize)

		% Bin QCFC data by distance and generate means and stds for each
		numThresholds = 10;
		BF_PlotQuantiles(allData(i).MeanEdgeWeight,allData(i).QCFCVec,numThresholds)
		hold on

		xlim([-1 2])
		ylim([-0.3 0.3])
		
		xlabel('Mean edge weight')
		ylabel('QC-FC')

		title(removeNoiseName,'Interpreter', 'none','FontSize',FSize,'FontWeight','normal')
	
		% ------------------------------------------------------------------------------
		% Plot: scatter of tSNR over mean FD
		% ------------------------------------------------------------------------------
		figure(Fig_QCFC_tSNR)
		subplot(4,ceil(numPrePro/4),i)
		hold on
		[tSNR,~] = GetTSNRForSample(allData(i).cfg,numSubs,numROIs,Parc);
		scatter(fdJenk_m,tSNR,'filled')

		% ylim([0 3000])
		ylim([0 1000])

		xlabel('Mean FD')
		ylabel('tSNR')
		
		title(removeNoiseName,'Interpreter', 'none','FontSize',FSize,'FontWeight','normal')
	end
end

% ------------------------------------------------------------------------------
% NBS contrasts
% ------------------------------------------------------------------------------
if runNBSPlots
	% Plot which contrasts?
	WhichContrasts = 'both'; % 'first' 'second' 'both'

	switch WhichContrasts
		case 'first'
			numContrasts = 1;
			contrasts = 1;
			offsets = 0;
		case 'second'
			numContrasts = 1;
			contrasts = 2;
			offsets = numPrePro;
		case 'both'
			numContrasts = 2;
			contrasts = [1,2];
			offsets = [0 numPrePro];
	end

	% ------------------------------------------------------------------------------
	% Figure size params
	% ------------------------------------------------------------------------------
	plotWidth = 10.5; plotHeight = plotWidth - 2.5;
	figWidth = numPrePro * plotWidth;
	figHeight = numContrasts * plotHeight;

	% Initialise figures
	Fig_T_Matrix = figure('color','w', 'units', 'centimeters', 'pos', [0 0 figWidth figHeight], 'name',['Fig_T_Matrix']); box('on'); movegui(Fig_T_Matrix,'center');
	Fig_T_DistDep = figure('color','w', 'units', 'centimeters', 'pos', [0 0 figWidth figHeight], 'name',['Fig_T_DistDep']); box('on'); movegui(Fig_T_DistDep,'center');
	Fig_Covar = figure('color','w', 'units', 'centimeters', 'pos', [0 0 figWidth figHeight], 'name',['Fig_Covar']); box('on'); movegui(Fig_Covar,'center');
	Fig_Var = figure('color','w', 'units', 'centimeters', 'pos', [0 0 figWidth figHeight], 'name',['Fig_Var']); box('on'); movegui(Fig_Var,'center');
	Fig_pVar = figure('color','w', 'units', 'centimeters', 'pos', [0 0 figWidth figHeight], 'name',['Fig_pVar']); box('on'); movegui(Fig_pVar,'center');

	for i = 1:numPrePro
	    removeNoise = allData(i).noiseOptions;
	    removeNoiseName = allData(i).noiseOptionsNames;

	    for j = 1:numContrasts
	    	con = contrasts(j);
	    	offset = offsets(j);

		    % ------------------------------------------------------------------------------
		    % Plot: test stat over distance
		    % ------------------------------------------------------------------------------
			figure(Fig_T_DistDep)
			subplot(numContrasts,numPrePro,i+offset)

			% Bin FD_FC data by distance and generate means and stds for each
			numThresholds = 10;
			BF_PlotQuantiles(ROIDistVec(allData(i).NaNFilter),allData(i).NBS_statVec{con}(allData(i).NaNFilter),numThresholds)
			ylim([-1.5 1.5])

			title(removeNoiseName,'Interpreter', 'none')
			hold on

			% ------------------------------------------------------------------------------
			% Plot: test stat as matrix
			% ------------------------------------------------------------------------------
			figure(Fig_T_Matrix)
			subplot(numContrasts,numPrePro,i+offset)

			ax = gca;
			ax.FontSize = FSize;

			% get data for plot
			TMatrix = allData(i).NBS_statMat{con};
			SigMatrix = allData(i).NBS_sigMat{con};
			NaNFilterMat = allData(i).NaNFilterMat;

			% filter test stat to retain significant connections
			TMatrix(SigMatrix ~= 1) = 0;
			clear SigMatrix

			% reorganise data by ROI struct
			TMatrix = TMatrix(ROI_idx,ROI_idx);
			NaNFilterMat = NaNFilterMat(ROI_idx,ROI_idx);

			% create ROI filter
			NaNFilterROI = NaNFilterMat(:,1);

			% filter data
			TMatrix = TMatrix(NaNFilterMat);
			% compute number of ROIs after filtering
			numROIsFilt = sum(NaNFilterROI);
			% reshape back to square (filtering outputs vector)
			TMatrix = reshape(TMatrix,numROIsFilt,numROIsFilt);

			% filter sorted ROI structure
			ROIStruct_temp = ROIStruct_com(NaNFilterROI,:);

			TheTMatrixPlot(TMatrix,ROIStruct_temp)

			% title(removeNoiseName,'Interpreter', 'none')
			xlabel('Brain region')
			ylabel('Brain region')
			title(removeNoiseName,'Interpreter', 'none','FontSize',FSize+2,'FontWeight','normal')
			% title([removeNoiseName,' (',num2str(round(PropSig(i,con),2)),'%)'],'Interpreter', 'none')
			hold on
			
			% ------------------------------------------------------------------------------
			% Plot: covariance
			% ------------------------------------------------------------------------------
			figure(Fig_Covar)
			subplot(numContrasts,numPrePro,i+offset)

			data = cell(1,numGroups);
			for g = 1:numGroups
				% Get covariance for group g
				data{g} = allData(i).VarCovar(:,:,Group == g);
				% Mean over subjects
				data{g} = mean(data{g},3);
				% Flatten
				data{g} = LP_FlatMat(data{g});
				% Threshold
				data{g} = data{g}(allData(i).NBS_sigVec{j} == 1);
				% data{g} = data{g}(allData(i).NaNFilter);
			end
			BF_JitteredParallelScatter(data,1,1,0)
			set(gca,'XTick',[1,2],'XTickLabel',{'Control','Patients'})
			% xlabel(['Contrast ',num2str(j)])
			ylabel('Covariance of edges')
			title(removeNoiseName,'Interpreter', 'none')
			ylim([-10 20])

			% ------------------------------------------------------------------------------
			% Plot: time series variance
			% ------------------------------------------------------------------------------
			figure(Fig_Var)
			subplot(numContrasts,numPrePro,i+offset)
			NaNFilterMat = allData(i).NaNFilterMat;
			NaNFilterROI = NaNFilterMat(:,1);

			data = cell(1,numGroups);
			for g = 1:numGroups
				% Get covariance for group g
				data{g} = allData(i).VarCovar(:,:,Group == g);
				% Mean over subjects
				data{g} = mean(data{g},3);
				% Get diagonal
				data{g} = data{g}(eye(numROIs) == 1);
				% Threshold
				data{g} = data{g}(NaNFilterROI);
			end
			BF_JitteredParallelScatter(data,1,1,0)
			set(gca,'XTick',[1,2],'XTickLabel',{'Control','Patients'})
			% xlabel(['Contrast ',num2str(j)])
			ylabel('Variance of edges')
			title(removeNoiseName,'Interpreter', 'none')
			ylim([0 40])

			% ------------------------------------------------------------------------------
			% Plot: pooled variance
			% ------------------------------------------------------------------------------
			figure(Fig_pVar)
			subplot(numContrasts,numPrePro,i+offset)

			data = cell(1,numGroups);
			for g = 1:numGroups
				% Get covariance for group g
				data{g} = allData(i).Var(:,:,Group == g);
				% Mean over subjects
				data{g} = mean(data{g},3);
				% Flatten
				data{g} = LP_FlatMat(data{g});
				% Threshold
				data{g} = data{g}(allData(i).NBS_sigVec{j} == 1);
				% data{g} = data{g}(allData(i).NaNFilter);
			end
			BF_JitteredParallelScatter(data,1,1,0)
			set(gca,'XTick',[1,2],'XTickLabel',{'Control','Patients'})
			% xlabel(['Contrast ',num2str(j)])
			ylabel('Sqrt(prod(var)) of edges')
			title(removeNoiseName,'Interpreter', 'none')
			ylim([0 25])
		end
	end
end

% ------------------------------------------------------------------------------
% T-stat distributions
% ------------------------------------------------------------------------------
if runTPlots
	% ------------------------------------------------------------------------------
	% Figure size params
	% ------------------------------------------------------------------------------
	plotWidth = 10.5; plotHeight = plotWidth - 2.5;
	figWidth = numPrePro * plotWidth;
	figHeight = plotHeight;
	
	% Initialise figures
	Fig_T_Hist = figure('color','w', 'units', 'centimeters', 'pos', [0 0 figWidth figHeight], 'name',['Fig_T_Hist']); box('on'); movegui(Fig_T_Hist,'center');
	
	for i = 1:numPrePro
	    removeNoise = allData(i).noiseOptions;
	    removeNoiseName = allData(i).noiseOptionsNames;

		% ------------------------------------------------------------------------------
		% Plot: test stat histogram
		% ------------------------------------------------------------------------------
		figure(Fig_T_Hist)
		subplot(1,numPrePro,i)

		% histogram(allData(i).NBS_statVec{1}(allData(i).NaNFilter),'FaceColor','k')
		[f,xi] = ksdensity(allData(i).NBS_statVec{1}(allData(i).NaNFilter));
		plot(xi,f,'--k','LineWidth',1.5)
		title(removeNoiseName,'Interpreter', 'none')
		hold on
		xlim([-10 10])

	end
end

% ------------------------------------------------------------------------------
% GSR T-stat plots for ICA-AROMA
% Note, this section has been coded to work specifically to look at ICA-AROMA+2P +/-GSR
% ------------------------------------------------------------------------------
if runGSRPlots
	% ------------------------------------------------------------------------------
	% Redefine noise pipeline variables to only look at those +/- GSR
	% ------------------------------------------------------------------------------
	% noiseOptions = {'sICA-AROMA+2P','sICA-AROMA+2P+GSR'};
	% noiseOptionsNames = {'ICA-AROMA+2Phys','ICA-AROMA+2Phys+GSR'};

	noiseOptions = {'24P+8P','24P+8P+4GSR';...
					'24P+aCC','24P+aCC+4GSR';...
					'sICA-AROMA+2P','sICA-AROMA+2P+GSR'};
	noiseOptionsNames = {'24HMP+8Phys','24HMP+8Phys+4GSR';...
						'24HMP+aCompCor','24HMP+aCompCor+4GSR';...
						'ICA-AROMA+2Phys','ICA-AROMA+2Phys+GSR'};

	% ------------------------------------------------------------------------------
	% Set contrasts
	% ------------------------------------------------------------------------------
	contrasts = [1,2;2,1;1,2];

	numPreProGSR = size(noiseOptions,1);

	% ------------------------------------------------------------------------------
	% Figure size params
	% ------------------------------------------------------------------------------
	plotWidth = 10.5; plotHeight = plotWidth - 2.5;
	figWidth = numPreProGSR * plotWidth;
	figHeight = plotHeight;
	
	% Initialise figures
	Fig_GSR = figure('color','w', 'units', 'centimeters', 'pos', [0 0 23 30], 'name',['Fig_GSR']); box('on'); movegui(Fig_GSR,'center');
	% Fig_GSRT_Hist = figure('color','w', 'units', 'centimeters', 'pos', [0 0 figWidth figHeight], 'name',['Fig_GSRT_Hist']); box('on'); movegui(Fig_GSRT_Hist,'center');
	% Fig_GSRT_GroupDiff = figure('color','w', 'units', 'centimeters', 'pos', [0 0 figWidth figHeight], 'name',['Fig_GSRT_GroupDiff']); box('on'); movegui(Fig_GSRT_GroupDiff,'center');
	% Fig_GSRT_Denom = figure('color','w', 'units', 'centimeters', 'pos', [0 0 figWidth figHeight], 'name',['Fig_GSRT_Denom']); box('on'); movegui(Fig_GSRT_Denom,'center');
	Fig_GSR_Covar = figure('color','w', 'units', 'centimeters', 'pos', [0 0 figWidth 2*figHeight], 'name',['Fig_GSR_Covar']); box('on'); movegui(Fig_GSR_Covar,'center');
	Fig_GSR_pVar = figure('color','w', 'units', 'centimeters', 'pos', [0 0 figWidth 2*figHeight], 'name',['Fig_GSR_pVar']); box('on'); movegui(Fig_GSR_pVar,'center');

	% containers
	Filter = cell(numPreProGSR,1);
	FCVec = cell(numPreProGSR,2);
	MeanFCVec = cell(numPreProGSR,2);
	T = cell(numPreProGSR,2);
	GroupDiff = cell(numPreProGSR,2);
	Denominator = cell(numPreProGSR,2);
	Covar = cell(numPreProGSR,2);
	pVar = cell(numPreProGSR,2);
	
	SpearRank = cell(numPreProGSR,1);

	for i = 1:numPreProGSR
	    removeNoiseName = noiseOptionsNames{i,1};

		clear idx
		% index of GSR- pipeline
		idx(1) = strmatch(noiseOptions{i,1},{allData.noiseOptions},'exact');
		% index of GSR+ pipeline
		idx(2) = strmatch(noiseOptions{i,2},{allData.noiseOptions},'exact');

		% ------------------------------------------------------------------------------
		% Edge filters
		% ------------------------------------------------------------------------------
		% NaNFilter
		% Filter{i} = allData(idx(1)).NaNFilter;
		% Common edges across GSR+/- for con.
		% Filter{i} = allData(idx(1)).NBS_sigVec{contrasts(i,1)} .* allData(idx(2)).NBS_sigVec{contrasts(i,1)};
		% Edges that were sig with GSR- and lost sig with GSR+. note order of subtraction is critical, check t matrix plots
		% Filter{i} = allData(idx(1)).NBS_sigVec{contrasts(i,1)} - allData(idx(2)).NBS_sigVec{contrasts(i,1)};
		% Edges that become sig in contrast 2 with GSR+ compare to GSR- contrast
		Filter{i} = allData(idx(2)).NBS_sigVec{contrasts(i,2)} - allData(idx(1)).NBS_sigVec{contrasts(i,1)};
		
		% Set -1 to 0
		Filter{i}(Filter{i} == -1) = 0;
		Filter{i} = logical(Filter{i});

		% ------------------------------------------------------------------------------
		% Loop over GSR+/-
		% ------------------------------------------------------------------------------
		for j = 1:2
			% Edges
			FCVec{i,j} = allData(idx(j)).FCVec(:,Filter{i});
			% Mean edges by group
			MeanFCVec{i,j}{1} = mean(FCVec{i,j}(Group == 1,:));
			MeanFCVec{i,j}{2} = mean(FCVec{i,j}(Group == 2,:));

			% T.
			T{i,j} = allData(idx(j)).NBS_statVec{contrasts(i,1)}(Filter{i});

			% Group diff and t-stat denominator
			[~,GroupDiff{i,j},Denominator{i,j}] = GetTStats(FCVec{i,j},Group,false);

			% Loop over groups to get Covar and pVar for group plots
			for g = 1:numGroups
				% Get covariance for group g
				Covar{i,j}{g} = allData(idx(j)).VarCovar(:,:,Group == g);
				% Mean over subjects
				Covar{i,j}{g} = mean(Covar{i,j}{g},3);
				% Flatten
				Covar{i,j}{g} = LP_FlatMat(Covar{i,j}{g});
				% Threshold
				Covar{i,j}{g} = Covar{i,j}{g}(Filter{i});

				% Get pooled variance for group g
				pVar{i,j}{g} = allData(idx(j)).Var(:,:,Group == g);
				% Mean over subjects
				pVar{i,j}{g} = mean(pVar{i,j}{g},3);
				% Flatten
				pVar{i,j}{g} = LP_FlatMat(pVar{i,j}{g});
				% Threshold
				pVar{i,j}{g} = pVar{i,j}{g}(Filter{i});
			end
		end

		% ------------------------------------------------------------------------------
		% Compute Phi coefficients between pipelines with and without GSR
		% ------------------------------------------------------------------------------
		% Check whether the connections introduced in the opposite direction 
		% with GSR+ overlap with the original connections from GSR-.
		% This indicates whether the same connections just get shift downward or
		% whether new connections are introduced with GSR+
		nonGSRSigVec = logical(allData(idx(1)).NBS_sigVec{contrasts(i,1)});
		GSRSigVec = logical(allData(idx(2)).NBS_sigVec{contrasts(i,2)});
		Phi = corr(nonGSRSigVec,GSRSigVec,'type','Pearson')

		% get subject-level spearman rank
		R = corr(FCVec{i,1}',FCVec{i,2}','type','Spearman');
		SpearRank{i} = R(eye(numSubs) == 1);
		mean(SpearRank{i}(Group == 1))
		mean(SpearRank{i}(Group == 2))

		% ------------------------------------------------------------------------------
		% Figures
		% ------------------------------------------------------------------------------

		% ------------------------------------------------------------------------------
		% Plot: test stat dist
		% ------------------------------------------------------------------------------
		figure(Fig_GSR)
		% figure(Fig_GSRT_Hist)
		subplot(3,numPreProGSR,i)
		hold on

		% Non GSR
		[f,xi] = ksdensity(T{i,1});
		plot(xi,f,'-k','LineWidth',1.5)

		% GSR
		[f,xi] = ksdensity(T{i,2});
		plot(xi,f,'-r','LineWidth',1.5)

		title(removeNoiseName,'Interpreter', 'none','FontSize',FSize+2,'FontWeight','normal')
		xlabel('Edge t-stats')
		xlim([-5 5])
		ylim([0 1.5])
		legend('GSR-','GSR+','Location','northoutside','Orientation','horizontal')
		
		% get spearman rank correlation
		Corr = corr(T{i,1},T{i,2},'type','Spearman');
		yLimits = ylim;
		xLimits = xlim;
		text(xLimits(1)+abs(xLimits(1))*.10,yLimits(2)-abs(yLimits(2))*.10,['Spearman: ',num2str(round(Corr,2))])

		% ------------------------------------------------------------------------------
		% Plot: group diff
		% ------------------------------------------------------------------------------
		% figure(Fig_GSRT_GroupDiff)
		subplot(3,numPreProGSR,i+numPreProGSR)
		hold on

		% Non GSR
		[f,xi] = ksdensity(GroupDiff{i,1});
		plot(xi,f,'-k','LineWidth',1.5)

		% GSR
		[f,xi] = ksdensity(GroupDiff{i,2});
		plot(xi,f,'-r','LineWidth',1.5)
		
		title(removeNoiseName,'Interpreter', 'none','FontSize',FSize+2,'FontWeight','normal')
		xlabel('Edge mean group difference')
		xlim([-0.5 0.5])
		% ylim([0 10])
		legend('GSR-','GSR+','Location','northoutside','Orientation','horizontal')

		% get spearman rank correlation
		Corr = corr(GroupDiff{i,1},GroupDiff{i,2},'type','Spearman');
		yLimits = ylim;
		xLimits = xlim;
		text(xLimits(1)+abs(xLimits(1))*.10,yLimits(2)-abs(yLimits(2))*.10,['Spearman: ',num2str(round(Corr,2))])

		% ------------------------------------------------------------------------------
		% Plot: t denom
		% ------------------------------------------------------------------------------
		% figure(Fig_GSRT_Denom)
		subplot(3,numPreProGSR,i+numPreProGSR*2)
		hold on

		% Non GSR
		[f,xi] = ksdensity(Denominator{i,1});
		plot(xi,f,'-k','LineWidth',1.5)

		% GSR
		[f,xi] = ksdensity(Denominator{i,2});
		plot(xi,f,'-r','LineWidth',1.5)
		
		title(removeNoiseName,'Interpreter', 'none','FontSize',FSize+2,'FontWeight','normal')
		xlabel('Edge t-stat denominator')
		xlim([0.02 0.09])
		ylim([0 100])
		legend('GSR-','GSR+','Location','northoutside','Orientation','horizontal')

		% get spearman rank correlation
		Corr = corr(Denominator{i,1},Denominator{i,2},'type','Spearman');
		yLimits = ylim;
		xLimits = xlim;
		text(xLimits(1)+abs(xLimits(1))*.10,yLimits(2)-abs(yLimits(2))*.10,['Spearman: ',num2str(round(Corr,2))])

		% ------------------------------------------------------------------------------
		% Plot: covariance
		% ------------------------------------------------------------------------------
		figure(Fig_GSR_Covar)
		for g = 1:numGroups
			if g == 1
				offset = 0;
			elseif g == 2
				offset = numPreProGSR;
			end

			subplot(numGroups,numPreProGSR,i+offset)
			hold on

			% Non GSR
			[f,xi] = ksdensity(Covar{i,1}{g});
			plot(xi,f,'-k','LineWidth',1.5)

			% GSR
			[f,xi] = ksdensity(Covar{i,2}{g});
			plot(xi,f,'-r','LineWidth',1.5)
			
			if g == 1
				title([removeNoiseName,': HCs'],'Interpreter', 'none','FontSize',FSize+2,'FontWeight','normal')
			elseif g == 2
				title([removeNoiseName,': Patients'],'Interpreter', 'none','FontSize',FSize+2,'FontWeight','normal')
			end
			xlabel('Edge ts cov')
			xlim([-5 10])
			ylim([0 2])
			legend('GSR-','GSR+','Location','northoutside','Orientation','horizontal')

			% get spearman rank correlation
			Corr = corr(Covar{i,1}{g},Covar{i,2}{g},'type','Spearman');
			yLimits = ylim;
			xLimits = xlim;
			text(xLimits(2)-abs(xLimits(2))*.50,yLimits(2)-abs(yLimits(2))*.10,['Spearman: ',num2str(round(Corr,2))])
		end

		% ------------------------------------------------------------------------------
		% Plot: pooled variance
		% ------------------------------------------------------------------------------
		figure(Fig_GSR_pVar)
		for g = 1:numGroups
			if g == 1
				offset = 0;
			elseif g == 2
				offset = numPreProGSR;
			end

			subplot(numGroups,numPreProGSR,i+offset)
			hold on

			% Non GSR
			[f,xi] = ksdensity(pVar{i,1}{g});
			plot(xi,f,'-k','LineWidth',1.5)

			% GSR
			[f,xi] = ksdensity(pVar{i,2}{g});
			plot(xi,f,'-r','LineWidth',1.5)
			
			if g == 1
				title([removeNoiseName,': HCs'],'Interpreter', 'none','FontSize',FSize+2,'FontWeight','normal')
			elseif g == 2
				title([removeNoiseName,': Patients'],'Interpreter', 'none','FontSize',FSize+2,'FontWeight','normal')
			end
			xlabel('Edge ts sqrt(prod(var))')
			xlim([0 20])
			ylim([0 1])
			legend('GSR-','GSR+','Location','northoutside','Orientation','horizontal')

			% get spearman rank correlation
			Corr = corr(pVar{i,1}{g},pVar{i,2}{g},'type','Spearman');
			yLimits = ylim;
			xLimits = xlim;
			text(xLimits(2)-abs(xLimits(2))*.25,yLimits(2)-abs(yLimits(2))*.10,['Spearman: ',num2str(round(Corr,2))])
		end
	end
end

% ------------------------------------------------------------------------------
% Overlap: sig edges across pipelines
% ------------------------------------------------------------------------------
if runOverlapPlots
	% ------------------------------------------------------------------------------
	% 1) Pairwise overlap between significant networks
	% ------------------------------------------------------------------------------
	WhichOverlap = 'jaccard'; % 'jaccard' 'phi' 
	f = figure('color','w', 'units', 'centimeters', 'pos', [0 0 23 12], 'name',['NBS Overlap']); box('on'); movegui(f,'center');
	for i = 1:numContrasts
		if i == 1
			data = cellfun(@(x) x(2),{allData(:).NBS_sigVec}); data = cell2mat(data);
		elseif i == 2
			data = cellfun(@(x) x(1),{allData(:).NBS_sigVec}); data = cell2mat(data);
		end

		% Retain only pipelines showing an effect
		pipeFilter = any(data);
		data = data(:,pipeFilter);
		data = full(data);

		% Filter NaNs
		NaNFilter = [allData(1).NaNFilter];
		data = data(NaNFilter,:); 

		switch WhichOverlap
			case 'jaccard'
				mat = [];
				for j = 1:size(data,2)
					for k = 1:size(data,2)
						x = data(:,j);
						y = data(:,k);
						% Num of signifcant edges in intersection
						Ci = sum(x + y == 2);
						% Num of signifcant edges in union
						Cu = sum(x + y > 0);

						mat(j,k) = Ci/Cu;
						% % express Ci as a fraction of the smaller NBS network
						% if sum(x) < sum(y)
						% 	mat(j,k) = Ci/sum(x);
						% elseif sum(y) < sum(x)
						% 	mat(j,k) = Ci/sum(y);
						% elseif sum(x) == sum(y)
						% 	mat(j,k) = 1;
						% end
					end
				end
			case 'phi'
				mat = corr(data,'type','Pearson');
		end

		% ------------------------------------------------------------------------------
		% Plot
		% ------------------------------------------------------------------------------
		subplot(1,2,i)
		imagesc(mat)
		axis square
		axis tight
		% colormap([flipud(BF_getcmap('blues',9,0));1,1,1;BF_getcmap('reds',9,0)])
		colormap(BF_getcmap('reds',9,0))
		caxis([0 1])
		colorbar
		ax = gca;
		ax.FontSize = FSize;
		x = [1:length(mat)];
		ax.XTick = x; ax.YTick = x;
		% ax.XTickLabel = '';
		ax.XTickLabel = noiseOptionsNames(:,pipeFilter);
		ax.XTickLabelRotation = 45;
		ax.YTickLabel = noiseOptionsNames(:,pipeFilter);
    	ax.TickLength = [0,0];

		% plot values in lower triangle of matrix
		for i = 1:size(mat,1)
			for j = i:size(mat,2)
				if i ~= j
					text(i,j,num2str(mat(i,j),'%0.1f'),'HorizontalAlignment','center',...
						'Color','k','FontSize',FSize,'FontWeight','normal');
				end
			end
		end

		if i == 1
			title('A) SCZ>HC','FontSize',FSize,'FontWeight','normal')
		elseif i == 2
			title('B) HC>SCZ','FontSize',FSize,'FontWeight','normal')
		end
	end

	% ------------------------------------------------------------------------------
	% 2) Figure 10
	% ------------------------------------------------------------------------------
	% Neuromarvl colours
	% BLUE | GREEN | YELLOW | RED
	% #0029ff | #00ff75 | #ffff00 |#ff0000

	cd('/Users/lindenmp/Dropbox/Work/Papers/PhD_Chapters/rfMRI_denoise/Figures/Fig10_Neuromarvl/neuromarvl')
	% node coords
	T = table(ROI_Coords(:,1),ROI_Coords(:,2),ROI_Coords(:,3),'VariableNames',{'x','y','z'});
	writetable(T,'coordinates.txt')
	clear T

	% node attributes
	T = table(ones(numROIs,1),ROIStructID,'VariableNames',{'dummy','ROIStructID'});
	writetable(T,'attributes.txt')

	% A)
    f = figure('color','w', 'units', 'centimeters', 'pos', [0 0 11.5 15], 'name',['24P+8P+4GSR_con2']); box('on'); movegui(f,'center');
    % subplot(1,2,1)
	% 24P+8P+4GSR
	SigMatrix = full(allData(6).NBS_sigMat{2});
	[out,outPC,outNorm] = plotClassifiedEdges(SigMatrix,ROIStructID,2,ROILabels);

	ConMatrix = full(allData(6).NBS_statMat{2});
	ConMatrix(SigMatrix == 0) = 0;
	dlmwrite('24P+8P+4GSR_con2.txt',ConMatrix)
	% http://immersive.erc.monash.edu.au/neuromarvl/?save=0c084e38-89e8-4708-8df7-5d0b9c381d9a_130.194.90.77
    % http://bit.ly/2sXxZiU

	% B)
    f = figure('color','w', 'units', 'centimeters', 'pos', [0 0 11.5 15], 'name',['ICA-AROMA+2P_con1']); box('on'); movegui(f,'center');
    % subplot(1,2,2)
	% ICA-AROMA+2P
	SigMatrix = full(allData(13).NBS_sigMat{1});
	[out,outPC,outNorm] = plotClassifiedEdges(SigMatrix,ROIStructID,2,ROILabels);

	% connectivity matrix
	ConMatrix = full(allData(13).NBS_statMat{1});
	ConMatrix(SigMatrix == 0) = 0;
	dlmwrite('ICA-AROMA+2P_con1.txt',ConMatrix)
	% http://immersive.erc.monash.edu.au/neuromarvl/?save=0ec770e2-6846-4725-ae05-a10d35706adf_130.194.90.77
	% http://bit.ly/2rIwrWx

	% C) 
    f = figure('color','w', 'units', 'centimeters', 'pos', [0 0 11.5 15], 'name',['6P+2P+GSR_con2']); box('on'); movegui(f,'center');
    % subplot(1,2,2)
	% 6P+2P+GSR_con2
	SigMatrix = full(allData(3).NBS_sigMat{2});
	[out,outPC,outNorm] = plotClassifiedEdges(SigMatrix,ROIStructID,2,ROILabels);

	ConMatrix = full(allData(3).NBS_statMat{2});
	ConMatrix(SigMatrix == 0) = 0;
	dlmwrite('6P+2P+GSR_con2.txt',ConMatrix)
	% http://immersive.erc.monash.edu.au/neuromarvl/?save=72338374-995e-4a93-b540-837511ff6608_130.194.90.77
	% http://bit.ly/2rYaua9
	
	% D)
    f = figure('color','w', 'units', 'centimeters', 'pos', [0 0 11.5 15], 'name',['6P+2P+GSR_con1']); box('on'); movegui(f,'center');
    % subplot(1,2,1)
	% 6P+2P+GSR_con1
	SigMatrix = full(allData(3).NBS_sigMat{1});
	[out,outPC,outNorm] = plotClassifiedEdges(SigMatrix,ROIStructID,2,ROILabels);

	ConMatrix = full(allData(3).NBS_statMat{1});
	ConMatrix(SigMatrix == 0) = 0;
	dlmwrite('6P+2P+GSR_con1.txt',ConMatrix)
	% http://immersive.erc.monash.edu.au/neuromarvl/?save=1c100289-5d8b-45fb-bae6-b8db4dbf37a4_130.194.90.77
	% http://bit.ly/2t1VTdr

end

clear pval* vec temp FC*
