% ------------------------------------------------------------------------------
% This script performs the primary analyses for
% Parkes, Fulcher, Yucel, Fornito. An evaluation of the efficacy, reliability, 
% and sensitivity of motion correction strategies for resting-state functional MRI
% 
% Copyright (C) 2017, Linden Parkes <lindenparkes@gmail.com>,
% ------------------------------------------------------------------------------
clear all; close all; clc

% ------------------------------------------------------------------------------
% Parent dirs
% ------------------------------------------------------------------------------
% MASSIVE
% parentdir = '/home/lindenmp/kg98_scratch/Linden/ResProjects/';
% ROIDir = '/projects/kg98/Linden/ROIs/';
% ExtHDD = [parentdir,'rfMRI_denoise/'];
% where the prepro scripts are
% funcdir = '/home/lindenmp/kg98/Linden/Scripts/rs-fMRI/func/';

% Mac (local)
% parentdir = '~/Dropbox/Work/ResProjects/';
% ROIDir = '~/Dropbox/Work/ROIs/';
% ExtHDD = '/Volumes/USB_D2/ResProjects/rfMRI_denoise/';
% where the prepro scripts are
% funcdir = '~/Dropbox/Work/git/rs-fMRI/func/';

% For those conducting reproducibility checks via figshare download
parentdir = '~/Desktop/figshare/';
ROIDir = '~/Desktop/ROIs/';
ExtHDD = '~/Desktop/figshare/rfMRI_denoise/';
% where the prepro scripts are
funcdir = '~/Desktop/rs-fMRI/func/';

addpath(funcdir)

% ------------------------------------------------------------------------------
% Set string switches
% ------------------------------------------------------------------------------
Projects = {'Beijing_Zang','UCLA','OCDPG','NYU_2','COBRE','COBRE_HCTSA','UCLA_HCTSA','NAMIC_HCTSA','ML_SNS','OCDPG_DCM'};
WhichProject = Projects{end};

WhichParc = 'Gordon'; % 'Gordon' 'Power'

if ismember('OCDPG',WhichProject,'rows') | ismember('UCLA',WhichProject,'rows') | ismember('COBRE',WhichProject,'rows') | ismember('OCDPG_DCM',WhichProject,'rows')
	% WhichSplit = 'Diagnostic'; % 'Motion' 'Diagnostic'
	WhichSplit = 'Motion'; % 'Motion' 'Diagnostic'
	% Note, this impacts whether only HCs are retained for the analyses or whether patients are included
	% If WhichSplit = 'Diagnostic' then the QC-FC plots will be different from what is shown in the manuscript because the patients will be included
	% whereas in the manuscript these benchmarks are calculated just using the HCs. Hence the default WhichSplit = 'Motion' will reproduce the primary
	% QC-FC figures accurately
elseif ismember('Beijing_Zang',WhichProject,'rows')
	WhichSplit = 'Motion'; % Only motion for Beijing
end

% ------------------------------------------------------------------------------
% Set logical switches
% ------------------------------------------------------------------------------
runPrimaryPlot = true;
runPlot = false;
runOverlapPlots = false;

fprintf(1, 'Running rfMRI QC.\n\tDataset: %s\n\tParcelation: %s\n',WhichProject,WhichParc);
if ismember('OCDPG',WhichProject,'rows') | ismember('UCLA',WhichProject,'rows') | ismember('COBRE',WhichProject,'rows')
	fprintf(1,'\tNBS: %s\n',WhichSplit);
end

% ------------------------------------------------------------------------------
% Set project variables
% ------------------------------------------------------------------------------
switch WhichProject
	case 'OCDPG_DCM'
		projdir = [parentdir,'rfMRI_DCM/OCDPG/'];
		sublist = [projdir,'OCDPGe.csv'];
		datadir = [projdir,'data/'];
        preprostr = '/func/prepro/';

		TR = 2.5;
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

% Only do this if there is more than a single group
if numel(unique(metadata.Diagnosis)) > 1
	% Retain only group 1 (assumed to be HCs) and group 2 (assumed to be patients)
	% I do this because the OCDPG dataset has some PGs in it that need to be removed
	metadata(metadata.Diagnosis == 3,:) = [];
end

numGroups = numel(unique(metadata.Diagnosis));

% ------------------------------------------------------------------------------
% Exclusion and censoring stuff
% ------------------------------------------------------------------------------
mname = 'rp*txt';
[metadata.exclude,metadata.mov,metadata.fdJenk,metadata.fdJenk_m,metadata.fdPower,metadata.fdPower_m,metadata.dvars,metadata.JP12ScrubMask,metadata.JP14ScrubMask] = GetExcludeForSample(datadir,metadata.ParticipantID,TR,preprostr,mname);
fprintf(1, 'done\n');

% compute number of volumes using the length of fdJenk
% note, this is assumed to be same for all subjects!
numVols = length(metadata.fdJenk{1});

if ismember('OCDPG_DCM',WhichProject,'rows')
	% exclude an outlier identified by Jeg during generation of the phenotypes
    idx = strmatch('sub-100',metadata.ParticipantID);
    metadata(idx,:) = [];

	% Also exclude sub-084 because of brain mask issues
	idx = strmatch('sub-084',metadata.ParticipantID);
	metadata(idx,:) = [];
end

% compute numsubs
numSubs = size(metadata,1);

% ------------------------------------------------------------------------------
% Create demographic table
% ------------------------------------------------------------------------------
Age = zeros(numGroups,1); Age_SD = zeros(numGroups,1); Gender = zeros(numGroups,1); Gender_percent = zeros(numGroups,1); 
IQ = zeros(numGroups,1); IQ_SD = zeros(numGroups,1); Medicated = zeros(numGroups,1); Medicated_percent = zeros(numGroups,1); 
fdJenk_m = zeros(numGroups,1); fdJenk_m_SD = zeros(numGroups,1); 
for i = 1:numGroups
    logi = metadata.Diagnosis == i;
    Age(i) = mean(metadata(logi,:).Age); Age_SD(i) = std(metadata(logi,:).Age);
    Gender(i) = sum(metadata(logi,:).Gender == 1); Gender_percent(i) = Gender(i) / sum(logi);
    IQ(i) = mean(metadata(logi,:).IQ); IQ_SD(i) = std(metadata(logi,:).IQ);
    Medicated(i) = sum(metadata(logi,:).Medicated == 1); Medicated_percent(i) = Medicated(i) / sum(logi);
    fdJenk_m(i) = mean(metadata(logi,:).fdJenk_m); fdJenk_m_SD(i) = std(metadata(logi,:).fdJenk_m);
end

demo_table = table(Age,Age_SD,Gender,Gender_percent,IQ,IQ_SD,Medicated,Medicated_percent,fdJenk_m,fdJenk_m_SD)
clear Age* Gender* IQ* Medicated* fdJenk_m* 

% ------------------------------------------------------------------------------
% Plot Movement params
% ------------------------------------------------------------------------------
fd = metadata.fdJenk_m;
fprintf(1, 'Mean FD for sample w/ lenient: %f (%f) \n', round(mean(fd(~metadata.exclude(:,1))),2), round(std(fd(~metadata.exclude(:,1))),2))
fprintf(1, 'Mean FD for sample w/ stringent: %f (%f) \n', round(mean(fd(~metadata.exclude(:,2))),2), round(std(fd(~metadata.exclude(:,2))),2))
fprintf(1, 'Mean FD for sample w/ spike reg: %f (%f) \n', round(mean(fd(~metadata.exclude(:,3))),2), round(std(fd(~metadata.exclude(:,3))),2))
fprintf(1, 'Mean FD for sample w/ basic scrub: %f (%f) \n', round(mean(fd(~metadata.exclude(:,4))),2), round(std(fd(~metadata.exclude(:,4))),2))
fprintf(1, 'Mean FD for sample w/ optimized scrub: %f (%f) \n\n', round(mean(fd(~metadata.exclude(:,5))),2), round(std(fd(~metadata.exclude(:,5))),2))

if numGroups > 1
	fd_HC = fd(metadata.Diagnosis == 1);
	fprintf(1, 'Mean FD for HCs w/ lenient: %f (%f) \n', round(mean(fd_HC(~metadata.exclude(metadata.Diagnosis == 1,1))),2), round(std(fd_HC(~metadata.exclude(metadata.Diagnosis == 1,1))),2))
	fprintf(1, 'Mean FD for HCs w/ stringent: %f (%f) \n', round(mean(fd_HC(~metadata.exclude(metadata.Diagnosis == 1,2))),2), round(std(fd_HC(~metadata.exclude(metadata.Diagnosis == 1,2))),2))
	fprintf(1, 'Mean FD for HCs w/ spike reg: %f (%f) \n', round(mean(fd_HC(~metadata.exclude(metadata.Diagnosis == 1,3))),2), round(std(fd_HC(~metadata.exclude(metadata.Diagnosis == 1,3))),2))
	fprintf(1, 'Mean FD for HCs w/ basic scrub: %f (%f) \n', round(mean(fd_HC(~metadata.exclude(metadata.Diagnosis == 1,4))),2), round(std(fd_HC(~metadata.exclude(metadata.Diagnosis == 1,4))),2))
	fprintf(1, 'Mean FD for HCs w/ optimized scrub: %f (%f) \n\n', round(mean(fd_HC(~metadata.exclude(metadata.Diagnosis == 1,5))),2), round(std(fd_HC(~metadata.exclude(metadata.Diagnosis == 1,5))),2))

	
	fd_Pat = fd(metadata.Diagnosis == 2);
	fprintf(1, 'Mean FD for patients w/ lenient: %f (%f) \n', round(mean(fd_Pat(~metadata.exclude(metadata.Diagnosis == 2,1))),2), round(std(fd_Pat(~metadata.exclude(metadata.Diagnosis == 2,1))),2))
	fprintf(1, 'Mean FD for patients w/ stringent: %f (%f) \n', round(mean(fd_Pat(~metadata.exclude(metadata.Diagnosis == 2,2))),2), round(std(fd_Pat(~metadata.exclude(metadata.Diagnosis == 2,2))),2))
	fprintf(1, 'Mean FD for patients w/ spike reg: %f (%f) \n', round(mean(fd_Pat(~metadata.exclude(metadata.Diagnosis == 2,3))),2), round(std(fd_Pat(~metadata.exclude(metadata.Diagnosis == 2,3))),2))
	fprintf(1, 'Mean FD for patients w/ basic scrub: %f (%f) \n', round(mean(fd_Pat(~metadata.exclude(metadata.Diagnosis == 2,4))),2), round(std(fd_Pat(~metadata.exclude(metadata.Diagnosis == 2,4))),2))
	fprintf(1, 'Mean FD for patients w/ optimized scrub: %f (%f) \n\n', round(mean(fd_Pat(~metadata.exclude(metadata.Diagnosis == 2,5))),2), round(std(fd_Pat(~metadata.exclude(metadata.Diagnosis == 2,5))),2))


	[h,p,~,stats] = ttest2(fd_HC,fd_Pat,'Vartype','unequal');
	if p < 0.05
		fprintf(1, 'There is a significant group difference in mean FD. t-value = %s. p-value = %s\n',num2str(stats.tstat),num2str(p));
	end

	BF_JitteredParallelScatter({fd_HC,fd_Pat});
	ax = gca;
	ax.XTick = [1,2];
	ax.XTickLabel = {'Controls','Patients'};
	xlabel('Group')
	ylabel('Mean FD')
end

% ------------------------------------------------------------------------------
% Retain only HCs for remainder of analysis
% ------------------------------------------------------------------------------
if ismember('Motion',WhichSplit,'rows')
	metadata = metadata(metadata.Diagnosis == 1,:);
end

% ------------------------------------------------------------------------------
% Calculate ICC: FD
% ------------------------------------------------------------------------------
if ismember('NYU_2',WhichProject,'rows')
	fprintf(1, 'Performing ICC on mean FD\n');
	[metadata.exclude_wrt,~,~,metadata.fdJenk_m_wrt,~,~,~,~,~] = GetExcludeForSample(datadir,metadata.ParticipantID,TR,'/session_1/func_2/prepro/');
	[metadata.exclude_brt,~,~,metadata.fdJenk_m_brt,~,~,~,~,~] = GetExcludeForSample(datadir,metadata.ParticipantID,TR,'/session_2/func_1/prepro/');

	% initialise
	FD_ICCw = zeros(1,numConnections);
	FD_ICCb = zeros(1,numConnections);

	% 1) within session FD_ICC
	x = [metadata.fdJenk_m, metadata.fdJenk_m_wrt];
	FD_ICCw = GetICC(x);

	% 2) between session FD_ICC
	y = [metadata.fdJenk_m, metadata.fdJenk_m_brt];
	FD_ICCb = GetICC(y);

	% We also need to create special censoring-based exclusion that span all scan sessions
	temp = metadata.exclude + metadata.exclude_wrt + metadata.exclude_brt;
	temp(temp > 0) = 1;
	metadata.exclude_all = temp;

	temp = metadata.exclude + metadata.exclude_wrt;
	temp(temp > 0) = 1;
	metadata.exclude_wrt = temp;

	temp = metadata.exclude + metadata.exclude_brt;
	temp(temp > 0) = 1;
	metadata.exclude_brt = temp;

	fprintf(1, '...done \n');
end

% backup metadata
metadata_bak = metadata;

% ------------------------------------------------------------------------------
% Variables
% ------------------------------------------------------------------------------
	allData = struct('noiseOptions',noiseOptions,...
					'noiseOptionsNames',noiseOptionsNames,...
					'PercExcluded',[],...
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
					'QCFC_AbsMed',[],...
					'QCFC_DistDep',[],...
					'QCFC_DistDep_Pval',[],...
					'MeanEdgeWeight',[],...
					'tDOF',[],...
					'tDOF_mean',[],...
					'tDOF_std',[],...
					'tDOF_gmean',[],...
					'tDOF_gstd',[]);

	if ismember('OCDPG',WhichProject,'rows') | ismember('UCLA',WhichProject,'rows') | ismember('COBRE',WhichProject,'rows')
		allData(1).NBS_statMat = [];
		allData(1).NBS_statVec = [];

		allData(1).NBS_sigMat = [];
		allData(1).NBS_sigVec = [];

		allData(1).NaNFilterMat = [];

		allData(1).NBS_PropSig = [];
		
		allData(1).u05_PropSig = [];
		allData(1).u001_PropSig = [];
		allData(1).FDR_PropSig = [];
	end

	% Add test-retest variables if WhichProject == 'NYU_2'
	if ismember('NYU_2',WhichProject,'rows')
		allData(1).PercExcludedwrt = [];
		allData(1).PercExcludedbrt = [];
		allData(1).FCblVec = [];
		allData(1).FCwVec = [];
		allData(1).FCbVec = [];
		allData(1).ICCw = [];
		allData(1).ICCw_mean = [];
		allData(1).ICCw_std = [];
		allData(1).ICCw_AbsMed = [];
		allData(1).ICCb = [];
		allData(1).ICCb_mean = [];
		allData(1).ICCb_std = [];
		allData(1).ICCb_AbsMed = [];
		allData(1).ICC_tstat = [];
		allData(1).ICC_Pval = [];
		allData(1).ICC_Pval_corr = [];
	end

% ------------------------------------------------------------------------------
% 						Loop over preprocessing pipelines
% ------------------------------------------------------------------------------
for i = 1:numPrePro
    WhichNoise = allData(i).noiseOptions;
    WhichNoiseName = allData(i).noiseOptionsNames;
	fprintf(1, '\nProcessing data: %s\n',WhichNoise);

    WhichNoiseSplit = strsplit(WhichNoise,'+');

	% ------------------------------------------------------------------------------
	% Get time series and functional connectivity data
	% ------------------------------------------------------------------------------
	cfgFile = 'cfg.mat';
	% WhichExclude = 5;
	if any(strmatch('JP12Scrub',WhichNoiseSplit,'exact')) == 1
		WhichNoise_temp = WhichNoiseSplit;
		WhichNoise_temp(end) = [];
		WhichNoise_temp = strjoin(WhichNoise_temp,'+');

		WhichExclude = 4;
		metadata = metadata_bak(~metadata_bak.exclude(:,WhichExclude),:);
		numSubs = size(metadata,1);
		fprintf(1, 'Excluded %u subjects \n', sum(metadata_bak.exclude(:,WhichExclude)));
		[allData(i).cfg,allData(i).FC,allData(i).FCVec,allData(i).VarCovar,allData(i).Var,allData(i).GCOR] = GetFCForSample(datadir,metadata.ParticipantID,preprostr,WhichNoise_temp,cfgFile,Parc,numROIs,numConnections,metadata.JP12ScrubMask);
	elseif any(strmatch('JP14Scrub',WhichNoiseSplit,'exact')) == 1
		WhichExclude = 5;
		metadata = metadata_bak(~metadata_bak.exclude(:,WhichExclude),:);
		numSubs = size(metadata,1);
		fprintf(1, 'Excluded %u subjects \n', sum(metadata_bak.exclude(:,WhichExclude)));
		[allData(i).cfg,allData(i).FC,allData(i).FCVec,allData(i).VarCovar,allData(i).Var,allData(i).GCOR] = GetFCForSample(datadir,metadata.ParticipantID,preprostr,WhichNoise,cfgFile,Parc,numROIs,numConnections,metadata.JP14ScrubMask);
	elseif any(strmatch('SpikeReg',WhichNoiseSplit,'exact')) == 1
		WhichExclude = 3;
		metadata = metadata_bak(~metadata_bak.exclude(:,WhichExclude),:);
		numSubs = size(metadata,1);
		fprintf(1, 'Excluded %u subjects \n', sum(metadata_bak.exclude(:,WhichExclude)));
		[allData(i).cfg,allData(i).FC,allData(i).FCVec,allData(i).VarCovar,allData(i).Var,allData(i).GCOR] = GetFCForSample(datadir,metadata.ParticipantID,preprostr,WhichNoise,cfgFile,Parc,numROIs,numConnections);
	else
		% For liberal
		WhichExclude = 1;
		% For stringent
		% WhichExclude = 2;
		metadata = metadata_bak(~metadata_bak.exclude(:,WhichExclude),:);
		numSubs = size(metadata,1);
		fprintf(1, 'Excluded %u subjects \n', sum(metadata_bak.exclude(:,WhichExclude)));
		[allData(i).cfg,allData(i).FC,allData(i).FCVec,allData(i).VarCovar,allData(i).Var,allData(i).GCOR] = GetFCForSample(datadir,metadata.ParticipantID,preprostr,WhichNoise,cfgFile,Parc,numROIs,numConnections);
	end

	% Store percentage of participants excluded
	allData(i).PercExcluded = sum(metadata_bak.exclude(:,WhichExclude) / size(metadata_bak,1)) * 100;

	% ------------------------------------------------------------------------------
	% Compute QC-FC
	% Note, make sure use z transformed r values
	% ------------------------------------------------------------------------------
	fprintf(1, 'Computing QC-FC: %s\n',WhichNoise);
	[allData(i).QCFCVec,allData(i).NaNFilter,allData(i).QCFC_PropSig_corr,allData(i).QCFC_PropSig_unc,allData(i).QCFC_AbsMed,allData(i).QCFC_DistDep,allData(i).QCFC_DistDep_Pval] = RunQCFC(metadata.fdJenk_m,allData(i).FC,ROIDistVec);

	% ------------------------------------------------------------------------------
	% Compute mean edge weight
	% ------------------------------------------------------------------------------
	fprintf(1, 'Computing mean edge weights: %s\n',WhichNoise);
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
	% Get tDOF
	% ------------------------------------------------------------------------------
	fprintf(1, 'Computing tDOF: %s\n',WhichNoise);
	allData(i).tDOF = zeros(numSubs,1);
	for j = 1:numSubs

		% get tDOF
		% First, find size of second dimension of noiseTS
		allData(i).tDOF(j) = size(allData(i).cfg(j).noiseTS,2);

		if any(strmatch('JP12Scrub',WhichNoiseSplit,'exact')) == 1
			allData(i).tDOF(j) = allData(i).tDOF(j) + sum(metadata.JP12ScrubMask{j});
		elseif any(strmatch('JP14Scrub',WhichNoiseSplit,'exact')) == 1
			allData(i).tDOF(j) = allData(i).tDOF(j) + sum(metadata.JP14ScrubMask{j});
		end

		% Then, if ICA-AROMA pipeline, find number of ICs and add to tDOF
		if ~isempty(strfind(WhichNoise,'ICA-AROMA'))
			x = dlmread([datadir,metadata.ParticipantID{j},preprostr,'ICA-AROMA_output/classified_motion_ICs.txt']);
			allData(i).tDOF(j) = allData(i).tDOF(j) + length(x);
		end
	end

	% ------------------------------------------------------------------------------
	% Calculate mean temporal degrees of freedom lost
	% ------------------------------------------------------------------------------
	% tDOF will be the same for most pipelines
	% but some have variable regressor amounts
	% so we take mean over subjects
	allData(i).tDOF_mean = mean(allData(i).tDOF);
	allData(i).tDOF_std = std(allData(i).tDOF);

	% also get mean by diagnostic groups (metadata.Diagnosis)
	if ismember('OCDPG',WhichProject,'rows') | ismember('UCLA',WhichProject,'rows') | ismember('COBRE',WhichProject,'rows')
		allData(i).tDOF_gmean(1) = mean(allData(i).tDOF(metadata.Diagnosis == 1));
		allData(i).tDOF_gmean(2) = mean(allData(i).tDOF(metadata.Diagnosis == 2));
		
		allData(i).tDOF_gstd(1) = std(allData(i).tDOF(metadata.Diagnosis == 1));
		allData(i).tDOF_gstd(2) = std(allData(i).tDOF(metadata.Diagnosis == 2));
	end

	% ------------------------------------------------------------------------------
	% Perform t-test on tDOF-loss
	% ------------------------------------------------------------------------------
	if ismember('OCDPG',WhichProject,'rows') | ismember('UCLA',WhichProject,'rows') | ismember('COBRE',WhichProject,'rows')
		x = allData(i).tDOF(metadata.Diagnosis == 1);
		y = allData(i).tDOF(metadata.Diagnosis == 2);

		[h,p,~,stats] = ttest2(x,y,'Vartype','unequal');
		if p < 0.05
			fprintf(1, 'Significant group difference in tDOF-loss. t-value(%s) = %s. p-value = %s\n',num2str(stats.df),num2str(stats.tstat),num2str(p));
		else
			fprintf(1, 'NO significant group difference in tDOF-loss. t-value(%s) = %s. p-value = %s\n',num2str(stats.df),num2str(stats.tstat),num2str(p));
		end
		fprintf(1, '\tMean tDOF-loss, group 1: %s\n', num2str(round(allData(i).tDOF_gmean(1),2)));
		fprintf(1, '\tMean tDOF-loss, group 2: %s\n', num2str(round(allData(i).tDOF_gmean(2),2)));
	end

	% ------------------------------------------------------------------------------
	% Get NBS contrasts
	% ------------------------------------------------------------------------------
	if ismember('OCDPG',WhichProject,'rows') | ismember('UCLA',WhichProject,'rows') | ismember('Beijing_Zang',WhichProject,'rows')
		fprintf(1, 'Getting NBS contrasts: %s\n',WhichNoise);
	    WhichStat = 'F';
	    numContrasts = 2;

	    switch WhichStat
	    	case 'F'
				if WhichExclude == 1
					nbsdir = [ExtHDD,WhichProject,'/NBS_liberal/'];
				elseif WhichExclude == 2
					nbsdir = [ExtHDD,WhichProject,'/NBS_stringent/'];
				else
					nbsdir = [ExtHDD,WhichProject,'/NBS_liberal/'];
				end

				if ismember('Beijing_Zang',WhichProject,'rows')
					nbsdir = [ExtHDD,WhichProject,'/NBS_liberal/'];
				end

			    % Load NBS data
			    file = dir([nbsdir,WhichSplit,'_',WhichParc,'_',WhichNoise,'_NBS*.mat']);
			    load([nbsdir,file(1).name],'nbsOut')

				nbs = nbsOut{1};
				% 1) generate residuals from NBS GLM
				% Get GLM
				GLM = nbs.GLM;
				% Number of independent GLM's to fit
				M = size(GLM.y,2);
				% Number of observations
				n = size(GLM.y,1);
				% Determine nuisance predictors not in contrast
				ind_nuisance = find(~GLM.contrast);
				% Regress out nuisance predictors and compute residual
				b = zeros(length(ind_nuisance),M);
				resid_y = zeros(n,M); 
				b = GLM.X(:,ind_nuisance) \ GLM.y;
				resid_y = GLM.y - GLM.X(:,ind_nuisance) * b; 

				fprintf(1, '\tFound %u networks\n', length(nbs.NBS.con_mat));
			
				for j = 1:numContrasts
					% 2) Do t-test using residuals
					x = resid_y(GLM.X(:,2) == 1,:);
					y = resid_y(GLM.X(:,2) == 0,:);
					if j == 1
						% G1 > G2 (HC > Patients; Low motion > high motion)
						[~,~,~,stats] = ttest2(x,y);
					elseif j == 2
						% G1 < G2 (HC < Patients; Low motion < high motion)
						[~,~,~,stats] = ttest2(y,x);
					end
					% Convert t-stat to square matrix
					tMat = LP_SquareVec(stats.tstat,numROIs);
					tMat = tMat + tMat';

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

				    % matrix of test statistics
					allData(i).NBS_statMat{j} = tMat;

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
						allData(i).NBS_sigMat{j} = full(nbs.NBS.con_mat{1});
						% get bottom triangle for accurate reordering of rows and columns (NBS only outputs upper triangle)
						allData(i).NBS_sigMat{j} = [allData(i).NBS_sigMat{j} + allData(i).NBS_sigMat{j}'];

						% Filter the sigMat by positive connections
						s = sign(allData(i).NBS_statMat{j});
						allData(i).NBS_sigMat{j}(s == -1) = 0;

						% flatten upper triangle
						allData(i).NBS_sigVec{j} = LP_FlatMat(allData(i).NBS_sigMat{j});
					end

					% ------------------------------------------------------------------------------
					% Calculate proportion of significant connections adjusted by NaNFilter
					% ------------------------------------------------------------------------------
				    allData(i).NBS_PropSig(j) = sum(allData(i).NBS_sigVec{j}(allData(i).NaNFilter)) / sum(allData(i).NaNFilter) * 100;

					% ------------------------------------------------------------------------------
					% Calculate PropSig at different thresholds for comparison
					% ------------------------------------------------------------------------------
					% Get design matrix w/ Cov from NBS
				    X = nbs.GLM.X;

				    % number of subjects in NBS
				    numSubsNBS = size(X,1);

					% Calculate p-value
					pval = 1 - tcdf(allData(i).NBS_statVec{j},numSubsNBS - rank(X));

					% 1) p < 0.05 unc
					pval_thr = pval < (0.05/2);
				    allData(i).u05_PropSig(j) = sum(pval_thr(allData(i).NaNFilter)) / sum(allData(i).NaNFilter) * 100;

					% 2) p < 0.001 unc
					pval_thr = pval < (0.001/2);
				    allData(i).u001_PropSig(j) = sum(pval_thr(allData(i).NaNFilter)) / sum(allData(i).NaNFilter) * 100;

					% 3) p < 0.05 FDR
					pval_cor = mafdr(pval,'BHFDR','true');
					pval_thr = pval_cor < (0.05/2);
				    allData(i).FDR_PropSig(j) = sum(pval_thr(allData(i).NaNFilter)) / sum(allData(i).NaNFilter) * 100;
				end
	    	case 'T'
				if WhichExclude == 1
					nbsdir = [ExtHDD,WhichProject,'/NBS_output_liberal/'];
				elseif WhichExclude == 2
					nbsdir = [ExtHDD,WhichProject,'/NBS_output_stringent/'];
				else
					nbsdir = [ExtHDD,WhichProject,'/NBS_output_liberal/'];
				end

				if ismember('Beijing_Zang',WhichProject,'rows')
					nbsdir = [ExtHDD,WhichProject,'/NBS_output_liberal/'];
				end

			    % Load NBS data
			    file = dir([nbsdir,WhichSplit,'_',WhichParc,'_',WhichNoise,'_NBS*.mat']);
			    load([nbsdir,file(1).name],'nbsOut')

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

					% ------------------------------------------------------------------------------
					% Calculate PropSig at different thresholds for comparison
					% ------------------------------------------------------------------------------
					% Get design matrix w/ Cov from NBS
				    X = nbs.GLM.X;

				    % number of subjects in NBS
				    numSubsNBS = size(X,1);

					% Calculate p-value
					pval = 1 - tcdf(allData(i).NBS_statVec{j},numSubsNBS - rank(X));

					% 1) p < 0.05 unc
					pval_thr = pval < (0.05/2);
				    allData(i).u05_PropSig(j) = sum(pval_thr(allData(i).NaNFilter)) / sum(allData(i).NaNFilter) * 100;

					% 2) p < 0.001 unc
					pval_thr = pval < (0.001/2);
				    allData(i).u001_PropSig(j) = sum(pval_thr(allData(i).NaNFilter)) / sum(allData(i).NaNFilter) * 100;

					% 3) p < 0.05 FDR
					pval_cor = mafdr(pval,'BHFDR','true');
					pval_thr = pval_cor < (0.05/2);
				    allData(i).FDR_PropSig(j) = sum(pval_thr(allData(i).NaNFilter)) / sum(allData(i).NaNFilter) * 100;
				end
		end

		% convert PropSig to full double. (for some reason it comes sparse)
		% allData(i).NBS_PropSig = full(allData(i).NBS_PropSig);

		% clean up
		clear nbsOut nbs
	end

	% ------------------------------------------------------------------------------
	% Compute ICC
	% ------------------------------------------------------------------------------
	if ismember('NYU_2',WhichProject,'rows')
		fprintf(1, 'Performing ICC for: %s\n',WhichNoise);

		% ------------------------------------------------------------------------------
		% Within session retest data
		% Get time series and functional connectivity data
		% ------------------------------------------------------------------------------
		% exclude_temp = metadata_bak.exclude_wrt;
		exclude_temp = metadata_bak.exclude_all;

		if any(strmatch('JP12Scrub',WhichNoiseSplit,'exact')) == 1
			WhichExclude = 4;
			metadata = metadata_bak(~exclude_temp(:,WhichExclude),:);
			fprintf(1, 'Excluded %u subjects \n', sum(exclude_temp(:,WhichExclude)));
			% Baseline under wrt exclusion
			[~,~,allData(i).FCblVec,~,~] = GetFCForSample(datadir,metadata.ParticipantID,'/session_1/func_1/prepro/',WhichNoise_temp,cfgFile,Parc,numROIs,numConnections,metadata.JP12ScrubMask);
			% wrt
			[~,~,allData(i).FCwVec,~,~] = GetFCForSample(datadir,metadata.ParticipantID,'/session_1/func_2/prepro/',WhichNoise_temp,cfgFile,Parc,numROIs,numConnections,metadata.JP12ScrubMask);
		elseif any(strmatch('JP14Scrub',WhichNoiseSplit,'exact')) == 1
			WhichExclude = 5;
			metadata = metadata_bak(~exclude_temp(:,WhichExclude),:);
			fprintf(1, 'Excluded %u subjects \n', sum(exclude_temp(:,WhichExclude)));
			% Baseline under wrt exclusion
			[~,~,allData(i).FCblVec,~,~] = GetFCForSample(datadir,metadata.ParticipantID,'/session_1/func_1/prepro/',WhichNoise,cfgFile,Parc,numROIs,numConnections,metadata.JP14ScrubMask);
			% wrt
			[~,~,allData(i).FCwVec,~,~] = GetFCForSample(datadir,metadata.ParticipantID,'/session_1/func_2/prepro/',WhichNoise,cfgFile,Parc,numROIs,numConnections,metadata.JP14ScrubMask);
		elseif any(strmatch('SpikeReg',WhichNoiseSplit,'exact')) == 1
			WhichExclude = 3;
			metadata = metadata_bak(~exclude_temp(:,WhichExclude),:);
			fprintf(1, 'Excluded %u subjects \n', sum(exclude_temp(:,WhichExclude)));
			% Baseline under wrt exclusion
			[~,~,allData(i).FCblVec,~,~] = GetFCForSample(datadir,metadata.ParticipantID,'/session_1/func_1/prepro/',WhichNoise,cfgFile,Parc,numROIs,numConnections);
			% wrt
			[~,~,allData(i).FCwVec,~,~] = GetFCForSample(datadir,metadata.ParticipantID,'/session_1/func_2/prepro/',WhichNoise,cfgFile,Parc,numROIs,numConnections);
		else
			WhichExclude = 1;
			% WhichExclude = 2;
			metadata = metadata_bak(~exclude_temp(:,WhichExclude),:);
			fprintf(1, 'Excluded %u subjects \n', sum(exclude_temp(:,WhichExclude)));
			% Baseline under wrt exclusion
			[~,~,allData(i).FCblVec,~,~] = GetFCForSample(datadir,metadata.ParticipantID,'/session_1/func_1/prepro/',WhichNoise,cfgFile,Parc,numROIs,numConnections);
			% wrt
			[~,~,allData(i).FCwVec,~,~] = GetFCForSample(datadir,metadata.ParticipantID,'/session_1/func_2/prepro/',WhichNoise,cfgFile,Parc,numROIs,numConnections);
		end

		% Store percentage of participants excluded
		allData(i).PercExcludedwrt = sum(exclude_temp(:,WhichExclude) / size(metadata_bak,1)) * 100;

		% ------------------------------------------------------------------------------
		% Calculate ICC
		% ------------------------------------------------------------------------------
		% initialise
		allData(i).ICCw = zeros(1,numConnections);

		% For each functional connection
		for j = 1:numConnections
			% 1) within session ICC
			x = [allData(i).FCblVec(:,j), allData(i).FCwVec(:,j)];
			allData(i).ICCw(j) = GetICC(x);
		end

		% ------------------------------------------------------------------------------
		% Between session retest data
		% Get time series and functional connectivity data
		% ------------------------------------------------------------------------------
		% exclude_temp = metadata_bak.exclude_brt;
		exclude_temp = metadata_bak.exclude_all;

		if any(strmatch('JP12Scrub',WhichNoiseSplit,'exact')) == 1
			WhichExclude = 4;
			metadata = metadata_bak(~exclude_temp(:,WhichExclude),:);
			fprintf(1, 'Excluded %u subjects \n', sum(exclude_temp(:,WhichExclude)));
			% Baseline under brt exclusion
			[~,~,allData(i).FCblVec,~,~] = GetFCForSample(datadir,metadata.ParticipantID,'/session_1/func_1/prepro/',WhichNoise_temp,cfgFile,Parc,numROIs,numConnections,metadata.JP12ScrubMask);
			% brt
			[~,~,allData(i).FCbVec,~,~] = GetFCForSample(datadir,metadata.ParticipantID,'/session_2/func_1/prepro/',WhichNoise_temp,cfgFile,Parc,numROIs,numConnections,metadata.JP12ScrubMask);
		elseif any(strmatch('JP14Scrub',WhichNoiseSplit,'exact')) == 1
			WhichExclude = 5;
			metadata = metadata_bak(~exclude_temp(:,WhichExclude),:);
			fprintf(1, 'Excluded %u subjects \n', sum(exclude_temp(:,WhichExclude)));
			% Baseline under brt exclusion
			[~,~,allData(i).FCblVec,~,~] = GetFCForSample(datadir,metadata.ParticipantID,'/session_1/func_1/prepro/',WhichNoise,cfgFile,Parc,numROIs,numConnections,metadata.JP14ScrubMask);
			% brt
			[~,~,allData(i).FCbVec,~,~] = GetFCForSample(datadir,metadata.ParticipantID,'/session_2/func_1/prepro/',WhichNoise,cfgFile,Parc,numROIs,numConnections,metadata.JP14ScrubMask);
		elseif any(strmatch('SpikeReg',WhichNoiseSplit,'exact')) == 1
			WhichExclude = 3;
			metadata = metadata_bak(~exclude_temp(:,WhichExclude),:);
			fprintf(1, 'Excluded %u subjects \n', sum(exclude_temp(:,WhichExclude)));
			% Baseline under brt exclusion
			[~,~,allData(i).FCblVec,~,~] = GetFCForSample(datadir,metadata.ParticipantID,'/session_1/func_1/prepro/',WhichNoise,cfgFile,Parc,numROIs,numConnections);
			% brt
			[~,~,allData(i).FCbVec,~,~] = GetFCForSample(datadir,metadata.ParticipantID,'/session_2/func_1/prepro/',WhichNoise,cfgFile,Parc,numROIs,numConnections);
		else
			WhichExclude = 1;
			% WhichExclude = 2;
			metadata = metadata_bak(~exclude_temp(:,WhichExclude),:);
			fprintf(1, 'Excluded %u subjects \n', sum(exclude_temp(:,WhichExclude)));
			% Baseline under brt exclusion
			[~,~,allData(i).FCblVec,~,~] = GetFCForSample(datadir,metadata.ParticipantID,'/session_1/func_1/prepro/',WhichNoise,cfgFile,Parc,numROIs,numConnections);
			% brt
			[~,~,allData(i).FCbVec,~,~] = GetFCForSample(datadir,metadata.ParticipantID,'/session_2/func_1/prepro/',WhichNoise,cfgFile,Parc,numROIs,numConnections);
		end

		% Store percentage of participants excluded
		allData(i).PercExcludedbrt = sum(exclude_temp(:,WhichExclude) / size(metadata_bak,1)) * 100;

		% ------------------------------------------------------------------------------
		% Calculate ICC
		% ------------------------------------------------------------------------------
		% initialise
		allData(i).ICCb = zeros(1,numConnections);

		% For each functional connection
		for j = 1:numConnections
			% 2) between session ICC
			y = [allData(i).FCblVec(:,j), allData(i).FCbVec(:,j)];
			allData(i).ICCb(j) = GetICC(y);
		end

		% ------------------------------------------------------------------------------
		% Get summary stats
		% ------------------------------------------------------------------------------
		% filter nans
		allData(i).ICCw = allData(i).ICCw(allData(i).NaNFilter);
		allData(i).ICCb = allData(i).ICCb(allData(i).NaNFilter);

		% mean
		allData(i).ICCw_mean = nanmean(allData(i).ICCw);
		allData(i).ICCb_mean = nanmean(allData(i).ICCb);

		% std
		allData(i).ICCw_std = nanstd(allData(i).ICCw);
		allData(i).ICCb_std = nanstd(allData(i).ICCb);

		% median abs
		allData(i).ICCw_AbsMed = nanmedian(abs(allData(i).ICCw));
		allData(i).ICCb_AbsMed = nanmedian(abs(allData(i).ICCb));

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

	% ------------------------------------------------------------------------------
	% Get num components for aCompCor50
	% ------------------------------------------------------------------------------
	if any(strmatch('12P+aCC50',WhichNoise,'exact')) == 1
		for j = 1:numSubs
			aCC50_num_wm(j) = dlmread([datadir,metadata.ParticipantID{j},preprostr,'12P+aCC50/aCC_num_wm.txt']);
			aCC50_num_csf(j) = dlmread([datadir,metadata.ParticipantID{j},preprostr,'12P+aCC50/aCC_num_csf.txt']);
		end
		aCC50_num_wm_mean = mean(aCC50_num_wm);
		aCC50_num_wm_std = std(aCC50_num_wm);

		aCC50_num_csf_mean = mean(aCC50_num_csf);
		aCC50_num_csf_std = std(aCC50_num_csf);
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

clear extraParams
% ------------------------------------------------------------------------------
% Chart colors and line styles
% ------------------------------------------------------------------------------
% tempColors = num2cell([255,105,97;97,168,255;178,223,138;117,112,179]./255,2);  
tempColors = num2cell([255,105,97;97,168,255;178,223,138;117,112,179;255,179,71]./255,2);  
for i = 1:numPrePro
	strs = strsplit(allData(i).noiseOptions,'+');
    % colors
    if any(strmatch('6P',strs,'exact')) == 1; theColors{i} = tempColors{1}; end
    if any(strmatch('24P',strs,'exact')) == 1; theColors{i} = tempColors{2}; end
    if any(strmatch('aCC',strs,'exact')) == 1 || any(strmatch('aCC50',strs,'exact')) == 1; theColors{i} = tempColors{3}; end
    if any(strmatch('ICA-AROMA',strs,'exact')) == 1; theColors{i} = tempColors{4}; end
    
    if any(strmatch('SpikeReg',strs,'exact')) == 1 | ...
		any(strmatch('JP12Scrub',strs,'exact')) == 1 | ...
		any(strmatch('JP14Scrub',strs,'exact')) == 1
    	theColors{i} = tempColors{5};
    end
    
    % line style
    if any(strmatch('GSR',strs,'exact')) == 1 | any(strmatch('4GSR',strs,'exact')) == 1 | any(strmatch('2GSR',strs,'exact')) == 1
    	theLines{i} = ':';
    	% theLines{i} = '--';
    else
    	theLines{i} = '-';
    end
end

% Plot labels
x = {allData(:).noiseOptionsNames};
y = num2cell(100 - [allData.PercExcluded]);
xy = cell(x);
for i = 1:length(x)
	if y{i} ~= 100
		xy{i} = strcat(x{i}, ' (',num2str(y{i},'%0.0f'),'%)');
	elseif y{i} == 100
		xy{i} = x{i};
	end
end

% ------------------------------------------------------------------------------
% Primary QC-FC figures
% ------------------------------------------------------------------------------
if runPrimaryPlot
	% ------------------------------------------------------------------------------
	% QC-FC significant proportion
	% Corresponds to Figure 1 and Figure 4 in case of censoring
	% ------------------------------------------------------------------------------
		Fig_QCFC_Dist = figure('color','w', 'units', 'centimeters', 'pos', [0 0 16 9], 'name',['Fig_QCFC_Dist']); box('on'); movegui(Fig_QCFC_Dist,'center');
		sp = subplot(1,3,1);
		pos = get(sp,'Position');
		% set(gca,'Position',[pos(1)*2.75, pos(2)*1.2, pos(3)*1.4, pos(4)*1]); % [left bottom width height]
		set(gca,'Position',[pos(1)*3.25, pos(2)*1.2, pos(3)*1.2, pos(4)*1]); % [left bottom width height]

		% Create data
		% data = {[allData(:).QCFC_PropSig_corr]'};
		data = {[allData(:).QCFC_PropSig_unc]'};
		data_std = cell(1,length(data)); [data_std{:}] = deal(zeros(size(data{1})));

		% Create table
		T = table(data{1},'RowNames',{allData(:).noiseOptionsNames}','VariableNames',{'QCFC_PropSig'})

		% Create bar chart
		clear extraParams
		extraParams.xTickLabels = xy;
		extraParams.xLabel = ''; % 'Pipeline'
		% extraParams.yLabel = 'QC-FC (%)';
		extraParams.yLabel = 'QC-FC uncorrected (%)';
		extraParams.theColors = theColors;
		extraParams.theLines = theLines;
		extraParams.yLimits = [0 110];

		TheBarChart(data,data_std,false,extraParams)

		% ------------------------------------------------------------------------------
		% QCFC distributions
		% ------------------------------------------------------------------------------
		sp = subplot(1,3,3);
		pos = get(sp,'Position');
		set(gca,'Position',[pos(1)*1, pos(2)*1.2, pos(3)*1.4, pos(4)*1]); % [left bottom width height]
		data2 = {allData(:).QCFCVec};
		clear extraParams
		% extraParams.theLabels = {allData(:).noiseOptionsNames};
		extraParams.customSpot = '';
		extraParams.add0Line = true;
		extraParams.theColors = theColors;
		BF_JitteredParallelScatter(data2,1,1,0,extraParams);
		ax = gca;

		% Set axis stuff
		ax.FontSize = FSize;
		% ax.XTick = [1:size(data{1},1)];
		ax.XTick = [];
		ax.XTickLabel = [];
		ax.XLim = ([0 numPrePro+1]);
		if ismember('NYU_2',WhichProject,'rows') | ismember('OCDPG',WhichProject,'rows')
			ax.YLim = ([-0.8 1.25]);
		else
			ax.YLim = ([-0.6 1]);
		end
		ylabel('QC-FC (Pearson''s r)')

	    % add text
	    TextRotation = 0;
        strprec = '%0.2f';
		data3 = {allData(:).QCFC_AbsMed};

	    text(1:size(data3,2),repmat(ax.YLim(2) - ax.YLim(2)*.05,1,size(data3,2)),num2str([data3{1,:}]',strprec),... 
	    'HorizontalAlignment','right',... 
	    'VerticalAlignment','middle',...
	    'Color','black',...
	    'FontSize', FSize,...
	    'Rotation',TextRotation)

		view(90,90)

	% ------------------------------------------------------------------------------
	% QC-FC distance dependence
	% Corresponds to Figure 2 and Figure 4 in case of censoring
	% ------------------------------------------------------------------------------
	data = {[allData(:).QCFC_DistDep]'};

	data_std = cell(1,length(data)); [data_std{:}] = deal(zeros(size(data{1})));
	
	% Create table
	T = table(data{1},'RowNames',{allData(:).noiseOptionsNames}','VariableNames',{'QCFC_DistDep'})

	% Create bar chart
	extraParams.xTickLabels = xy;
	extraParams.xLabel = '';
	extraParams.yLabel = 'QC-FC distance dependence (Spearman''s rho)';
    extraParams.theColors = theColors;
    extraParams.theLines = theLines;
	extraParams.yLimits = [-0.5 0.5];

	Fig_QCFC_DistDep = figure('color','w', 'units', 'centimeters', 'pos', [0 0 10.5 9], 'name',['Fig_QCFC_DistDep']); box('on'); movegui(Fig_QCFC_DistDep,'center');
	sp = subplot(1,1,1);
	pos = get(sp,'Position');
	set(gca,'Position',[pos(1)*4.9, pos(2)*1.2, pos(3)*0.425, pos(4)*1]); % [left bottom width height]

	TheBarChart(data,data_std,false,extraParams)
end

% ------------------------------------------------------------------------------
% Figures
% ------------------------------------------------------------------------------
if runPlot
	% ------------------------------------------------------------------------------
	% QC-FC corrected
	% ------------------------------------------------------------------------------
	% Create data
	data = {[allData(:).QCFC_PropSig_corr]'};
	data_std = cell(1,length(data)); [data_std{:}] = deal(zeros(size(data{1})));

	% Create table
	T = table(data{1},'RowNames',{allData(:).noiseOptionsNames}','VariableNames',{'QCFC_PropSig_corr'})

	% Create bar chart
	clear extraParams
	% extraParams.xTickLabels = {allData(:).noiseOptionsNames};
	extraParams.xTickLabels = xy;
	extraParams.xLabel = ''; % 'Pipeline'
	extraParams.yLabel = 'QC-FC FDR-corr. (%)';
    extraParams.theColors = theColors;
    extraParams.theLines = theLines;
	extraParams.yLimits = [0 110];

	Fig_QCFC_FDR = figure('color','w', 'units', 'centimeters', 'pos', [0 0 10.5 9], 'name',['Fig_QCFC_FDR']); box('on'); movegui(Fig_QCFC_FDR,'center');
	sp = subplot(1,1,1);
	pos = get(sp,'Position');
	set(gca,'Position',[pos(1)*4.9, pos(2)*1.2, pos(3)*0.425, pos(4)*1]); % [left bottom width height]

	TheBarChart(data,data_std,false,extraParams)

	% ------------------------------------------------------------------------------
	% Figures: QCFC
	% ------------------------------------------------------------------------------
	% Initialise figures
	Fig_QCFC_DistDepBig = figure('color','w', 'units', 'centimeters', 'pos', [0 0 25 27], 'name',['Fig_QCFC_DistDepBig']); box('on'); movegui(Fig_QCFC_DistDepBig,'center');

	pipelines2Retain = logical([0 0 0 0 1 1 1 1 0 0 0 0 1 1 0 0 1 1 1]);
	noiseOptions_temp = {allData(pipelines2Retain).noiseOptions};
	noiseOptionsNames_temp = {allData(pipelines2Retain).noiseOptionsNames};
	theColors_temp = theColors(pipelines2Retain);

	numPrePro_temp = length(noiseOptions_temp);

	for i = 1:numPrePro_temp

	    WhichNoise = noiseOptions_temp{i};
	    WhichNoiseName = noiseOptionsNames_temp{i};
		idx = strmatch(WhichNoise,{allData(:).noiseOptions},'exact');

		% ------------------------------------------------------------------------------
		% Plot: distance dependence
		% ------------------------------------------------------------------------------
		figure(Fig_QCFC_DistDepBig)
		subplot(4,ceil(numPrePro_temp/4),i)
		set(gca,'FontSize',FSize)

		% Bin QCFC data by distance and generate means and stds for each
		numThresholds = 11;
		BF_PlotQuantiles(ROIDistVec(allData(idx).NaNFilter),allData(idx).QCFCVec,numThresholds,0,0,theColors_temp{i})
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

		title(WhichNoiseName,'Interpreter', 'none','FontSize',FSize,'FontWeight','normal')
	end

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
	extraParams.xTickLabels = xy;
	extraParams.xLabel = '';
	extraParams.yLabel = 'tDOF-loss (# regressors)';
    extraParams.theColors = theColors;
    extraParams.theLines = theLines;
	extraParams.yLimits = [0 140];

	Fig_tDOF = figure('color','w', 'units', 'centimeters', 'pos', [0 0 10.5 9], 'name',['Fig_tDOF']); box('on'); movegui(Fig_tDOF,'center');
	sp = subplot(1,1,1);
	pos = get(sp,'Position');
	set(gca,'Position',[pos(1)*4.9, pos(2)*1.2, pos(3)*0.425, pos(4)*1]); % [left bottom width height]

	TheBarChart(data,data_std,false,extraParams)

	% ------------------------------------------------------------------------------
	% Size of significant edge component
	% Corresponds to Figure 6
	% ------------------------------------------------------------------------------
	if ~ismember('NYU_2',WhichProject,'rows')
		% Create table
		if ismember('Diagnostic',WhichSplit,'rows')
			PropSig = vertcat(allData(:).NBS_PropSig);
			% PropSig = vertcat(allData(:).u05_PropSig);
			% PropSig = vertcat(allData(:).u001_PropSig);
			% PropSig = vertcat(allData(:).FDR_PropSig);
		elseif ismember('Motion',WhichSplit,'rows')
			% PropSig = vertcat(allData(:).NBS_PropSig);
			PropSig = vertcat(allData(:).u05_PropSig);
			% PropSig = vertcat(allData(:).u001_PropSig);
			% PropSig = vertcat(allData(:).FDR_PropSig);
		end

		data = {PropSig(:,1),PropSig(:,2)};
		data_std = cell(1,length(data)); [data_std{:}] = deal(zeros(size(data{1})));
		T = table([data{:}],'RowNames',{allData(:).noiseOptionsNames}','VariableNames',{'PropSig'})

		% Create bar chart
		extraParams.xTickLabels = xy;
		extraParams.xLabel = '';
		extraParams.yLabel = 'Significant connections (%)';
		extraParams.theColors = theColors;
		extraParams.theLines = theLines;
		extraParams.makeABS = true;
		if ismember('Diagnostic',WhichSplit,'rows')
			extraParams.yLimits = [-40 40];
		elseif ismember('Motion',WhichSplit,'rows')
			extraParams.yLimits = [-80 40];
			% extraParams.yLimits = [-30 10];
		end

		Fig_NBS = figure('color','w', 'units', 'centimeters', 'pos', [0 0 10.5 9], 'name',['Fig_NBS']); box('on'); movegui(Fig_NBS,'center');
		sp = subplot(1,1,1);
		pos = get(sp,'Position');
		set(gca,'Position',[pos(1)*4.9, pos(2)*1.2, pos(3)*0.44, pos(4)*1]); % [left bottom width height]

		TheBarChart(data,data_std,false,extraParams)
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
		sp = subplot(1,1,1);
		pos = get(sp,'Position');
		set(gca,'Position',[pos(1)*4.9, pos(2)*1.2, pos(3)*0.425, pos(4)*1]); % [left bottom width height]
		% set(gca,'Position',[pos(1)*4.65, pos(2)*1.2, pos(3)*0.475, pos(4)*1]); % [left bottom width height]

		data = {allData(:).ICCw}';

		% Plot labels
		x = {allData(:).noiseOptionsNames};
		y = num2cell(100 - [allData.PercExcludedwrt]);
		xy = cell(x);
		for i = 1:length(x)
			if y{i} ~= 100
				xy{i} = strcat(x{i}, ' (',num2str(y{i},'%0.0f'),'%)');
			elseif y{i} == 100
				xy{i} = x{i};
			end
		end

		clear extraParams
		extraParams.customSpot = '';
		extraParams.add0Line = true;
		extraParams.theColors = theColors;
		BF_JitteredParallelScatter(data,1,1,0,extraParams);
		ax = gca;

		% Set axis stuff
		ax.FontSize = FSize;
		ax.XTick = 1:numPrePro;
		ax.XTickLabel = xy;
		ax.XLim = ([0 numPrePro+1]);
		ax.YLim = ([-1 1.5]);
		xlabel('')
		ylabel('Within Session ICC')

	    % add text
	    TextRotation = 0;
        strprec = '%0.2f';
		data2 = {allData(:).ICCw_AbsMed};

	    text(1:size(data2,2),repmat(ax.YLim(2) - ax.YLim(2)*.05,1,size(data2,2)),num2str([data2{1,:}]',strprec),... 
	    'HorizontalAlignment','right',... 
	    'VerticalAlignment','middle',...
	    'Color','black',...
	    'FontSize', FSize,...
	    'Rotation',TextRotation)

		view(90,90)

		% ------------------------------------------------------------------------------
		% BRT ICC
		% ------------------------------------------------------------------------------
		Fig_ICCb = figure('color','w', 'units', 'centimeters', 'pos', [0 0 10.5 9], 'name',['Fig_ICCb']); box('on'); movegui(Fig_ICCb,'center');
		sp = subplot(1,1,1);
		pos = get(sp,'Position');
		set(gca,'Position',[pos(1)*4.9, pos(2)*1.2, pos(3)*0.425, pos(4)*1]); % [left bottom width height]
		% set(gca,'Position',[pos(1)*4.65, pos(2)*1.2, pos(3)*0.475, pos(4)*1]); % [left bottom width height]
		
		data = {allData(:).ICCb}';

		% Plot labels
		x = {allData(:).noiseOptionsNames};
		y = num2cell(100 - [allData.PercExcludedbrt]);
		xy = cell(x);
		for i = 1:length(x)
			if y{i} ~= 100
				xy{i} = strcat(x{i}, ' (',num2str(y{i},'%0.0f'),'%)');
			elseif y{i} == 100
				xy{i} = x{i};
			end
		end

		clear extraParams
		extraParams.customSpot = '';
		extraParams.add0Line = true;
		extraParams.theColors = theColors;
		BF_JitteredParallelScatter(data,1,1,0,extraParams);
		ax = gca;

		% Set axis stuff
		ax.FontSize = FSize;
		ax.XTick = 1:numPrePro;
		ax.XTickLabel = xy;
		ax.XLim = ([0 numPrePro+1]);
		ax.YLim = ([-1 1.5]);
		xlabel('')
		ylabel('Between Session ICC')

	    % add text
	    TextRotation = 0;
        strprec = '%0.2f';
		data2 = {allData(:).ICCb_AbsMed};

	    text(1:size(data2,2),repmat(ax.YLim(2) - ax.YLim(2)*.05,1,size(data2,2)),num2str([data2{1,:}]',strprec),... 
	    'HorizontalAlignment','right',... 
	    'VerticalAlignment','middle',...
	    'Color','black',...
	    'FontSize', FSize,...
	    'Rotation',TextRotation)

		view(90,90)
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
		ax.XTickLabel = xy(:,pipeFilter);
		ax.XTickLabelRotation = 45;
		ax.YTickLabel = xy(:,pipeFilter);
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
	% Figure 10: neuromarvl
	% Neuromarvl colours
	% BLUE | GREEN | YELLOW | RED
	% #0029ff | #00ff75 | #ffff00 | #ff0000
	% ------------------------------------------------------------------------------
	pipelines2Retain = logical([0 0 0 0 0 1 0 1 0 0 0 0 0 1 0 0 1 0 0]);
	noiseOptions_temp = {allData(pipelines2Retain).noiseOptions};

	% cd('/Users/lindenmp/Dropbox/Work/Papers/PhD_Chapters/rfMRI_denoise/Figures/Fig10_Neuromarvl/neuromarvl')
	% cd('/Users/lindenmp/Dropbox/Work/Papers/PhD_Chapters/rfMRI_denoise/Submissions/NeuroImage/Revision/Figures/Fig10_Neuromarvl')
	% node coords
	T = table(ROI_Coords(:,1),ROI_Coords(:,2),ROI_Coords(:,3),'VariableNames',{'x','y','z'});
	writetable(T,'coordinates.txt')
	clear T

	% node attributes
	T = table(ones(numROIs,1),ROIStructID,'VariableNames',{'dummy','ROIStructID'});
	writetable(T,'attributes.txt')

	% Neuromarvl files
	for i = 1:length(noiseOptions_temp)
		idx = strmatch(noiseOptions_temp{i},{allData(:).noiseOptions},'exact');
		for j = 1:numContrasts
			SigMatrix = full(allData(idx).NBS_sigMat{j});
			ConMatrix = full(allData(idx).NBS_statMat{j});
			ConMatrix(SigMatrix == 0) = 0;
			dlmwrite([noiseOptions_temp{i},'_con',num2str(j),'.txt'],ConMatrix)

			% Get t-value colourbar upper and lower (for illustrator)
			x = ConMatrix(:);
			x = sort(x,'descend');
			fprintf(1, '%s %s: T-value upper: %s lower: %s \n', noiseOptions_temp{i}, num2str(j), num2str(round(x(1),2)),num2str(round(x(200),2)));
		end
	end

	% ------------------------------------------------------------------------------
	% Figure 10: matrix plots
	% ------------------------------------------------------------------------------
	for i = 1:length(noiseOptions_temp)
		idx = strmatch(noiseOptions_temp{i},{allData(:).noiseOptions},'exact');
		
		% Contrast 1
		SigMatrix = full(allData(idx).NBS_sigMat{1});
		[~,plotMat,~] = plotClassifiedEdges(SigMatrix,ROIStructID);
		% add constant
		plotMat = plotMat + 1;
		% expand
		x = triu(plotMat);
		x = [zeros(numROIComms,1), x]; x = [x; zeros(1,numROIComms+1)];
		x = [zeros(numROIComms+1,1), x]; x = [x; zeros(1,numROIComms+2)];

		% Contrast 2
		SigMatrix = full(allData(idx).NBS_sigMat{2});
		[~,plotMat,~] = plotClassifiedEdges(SigMatrix,ROIStructID);
		% add constant
		plotMat = plotMat + 1;
		% expand
		y = tril(plotMat);
		y = [y, zeros(numROIComms,1)]; y = [zeros(1,numROIComms+1); y];
		y = [y, zeros(numROIComms+1,1)]; y = [zeros(1,numROIComms+2); y];

		% Combine
		plotMat = x + y;

    	% Create mask
	    imAlpha = ones(size(plotMat));
		imAlpha(plotMat == 0) = 0;
		
		% remove the constant
		plotMat = plotMat - 1;
		% convert to %
    	plotMat = plotMat * 100;

		f = figure('color','w', 'units', 'centimeters', 'pos', [0 0 7 7], 'name',[noiseOptions_temp{i},'_con',num2str(j)]); box('off');
	    cmax = 12; 
	    % cmax = ceil(max(plotMat(:)));
		imagesc(plotMat,'AlphaData',imAlpha);
	        axis square
	        axis tight
	        colormap(BF_getcmap('reds',6,0))
	        caxis([0 cmax])
		set(gca,'TickLength', [0 0])
		% set(gca,'XTick',''); % Set x tick vals
		set(gca,'XTick',1:1:numROIComms); % Set x tick vals
		set(gca,'XTickLabel',ROILabels,'FontSize',8,'FontWeight','normal','XTickLabelRotation',90); % set x tick labels. these will be off and need to be corrected in illustrator
		set(gca,'YTick',1:1:numROIComms); % set y tick vals
		set(gca,'YTickLabel',ROILabels,'FontSize',8,'FontWeight','normal'); % set y tick labels
		colorbar; % show colour bar

		set(f, 'PaperPositionMode', 'auto')
		print(f,noiseOptions_temp{i},'-dsvg')
	end
end

clear pval* vec temp FC*
