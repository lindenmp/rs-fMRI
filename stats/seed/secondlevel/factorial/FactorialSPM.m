% This script will run second level estimation for a resting state fMRI analysis.
% Note that this script assumes you have one factor (hemisphere) that has two levels (left and right)
% This script assumes that you have got paired con images that were output from 1st level analysis. 
% That is, one con image per hemisphere that correspond to identical ROIs from the left and right hemispheres
% These must be stored in separate folders (e.g., FirstLevel_L/ and FirstLevel_R/)
% 	e.g., <SUBJECTID>/rfMRI/FirstLevel_L/con_0001.nii & <SUBJECTID>/rfMRI/FirstLevel_R/con_0001.nii
%
% This script loops over subjects specificied in the corresponding input and takes all the con_000# images

% ------------------------------------------------------------------------------
% Linden Parkes, Brain & Mental Health Laboratory, 2017
% ------------------------------------------------------------------------------
clear all; close all; clc
rng('default')

% ------------------------------------------------------------------------------
% Parent dir
% ------------------------------------------------------------------------------
parentdir = '/home/lindenmp/kg98/Linden/';
parentdir_scratch = '/home/lindenmp/kg98_scratch/Linden/';

% ------------------------------------------------------------------------------
% Add paths - edit this section
% ------------------------------------------------------------------------------
% where the prepro scripts are
funcdir = [parentdir,'Scripts/rs-fMRI/func/'];
addpath(funcdir)

% where spm12 is
spmdir = [parentdir,'Scripts/Tools/spm12/'];
addpath(spmdir)

% ------------------------------------------------------------------------------
% Set options
% ------------------------------------------------------------------------------
Projects = {'OCDPG_DCM','GenCog','UCLA'};
WhichProject = Projects{1}

WhichSeed = 'TriStri'
% WhichSeed = 'DiMartino'
% WhichSeed = 'GSR'

% ------------------------------------------------------------------------------
% Set parcellation
% Note, code is not setup to process multiple parcellations concurrently.
% ------------------------------------------------------------------------------
switch WhichSeed
    case 'TriStri'
    	conNums = [2,3]; % 2 = dorsal, 3 = ventral
    case 'DiMartino'
    	conNums = [3,1]; % 3 = DC, 1 = VSi
    case 'GSR'
    	conNums = [1];
end

% ------------------------------------------------------------------------------
% Set project variables
% ------------------------------------------------------------------------------
switch WhichProject
    case 'OCDPG_DCM'
        projdir = [parentdir_scratch,'ResProjects/rfMRI_DCM/OCDPG/'];
        datadir = [projdir,'data/'];
        preprostr = '/func/prepro/';

        datafile = [projdir,'OCDPGe.csv'];
        TR = 2.5;

		% WhichNoise = '6P/'
		WhichNoise = 'ICA-AROMA+2P/'
		% WhichNoise = 'ICA-AROMA+2P+GSR/'
    case 'GenCog'
        projdir = [parentdir_scratch,'ResProjects/rfMRI_DCM/GenCog/'];
        datadir = [projdir,'data/'];
        preprostr = '/func/prepro/';

        datafile = [projdir,'GenCog.csv'];
        TR = 0.754;

		WhichNoise = 'ICA-FIX/'
		% WhichNoise = 'ICA-FIX+GSR/'
    case 'UCLA'
        projdir = [parentdir_scratch,'ResProjects/rfMRI_denoise/UCLA/'];
        datadir = [projdir,'data/'];
        preprostr = '/func/prepro/';

        datafile = [projdir,'UCLA.csv'];
        TR = 2;

		% WhichNoise = '6P/'
		WhichNoise = 'ICA-AROMA+2P/'
		% WhichNoise = 'ICA-AROMA+2P+GSR/'
end

% ------------------------------------------------------------------------------
% Non-imaging metadata
% ------------------------------------------------------------------------------
metadata = readtable(datafile);

if ismember('OCDPG_DCM',WhichProject,'rows')
	% exclude an outlier identified by Jeg during generation of the phenotypes
    idx = strmatch('sub-100',metadata.ParticipantID);
    metadata(idx,:) = [];

	% Also exclude sub-084 because of brain mask issues
	idx = strmatch('sub-084',metadata.ParticipantID);
	metadata(idx,:) = [];
end

numGroups = length(unique(metadata.Diagnosis));

% ------------------------------------------------------------------------------
% Exclusion and censoring stuff
% ------------------------------------------------------------------------------
if ismember('OCDPG_DCM',WhichProject,'rows') | ismember('UCLA',WhichProject,'rows')
	[metadata.exclude,metadata.mov,metadata.fdJenk,metadata.fdJenk_m,metadata.fdPower,metadata.fdPower_m,metadata.dvars,metadata.JP12ScrubMask,metadata.JP14ScrubMask] = GetExcludeForSample(datadir,metadata.ParticipantID,TR,preprostr);
elseif ismember('GenCog',WhichProject,'rows')
	[metadata.exclude,metadata.mov,metadata.fdJenk,metadata.fdJenk_m,metadata.fdPower,metadata.fdPower_m,metadata.dvars,metadata.JP12ScrubMask,metadata.JP14ScrubMask] = GetExcludeForSample(datadir,metadata.ParticipantID,TR,preprostr,'prefiltered_func_data_mcf.par');
end

fprintf(1, 'done\n');

% WhichExclude = 1;
WhichExclude = 2;
fprintf(1, 'Excluded %u subjects \n', sum(metadata.exclude(:,WhichExclude)));
metadata = metadata(~metadata.exclude(:,WhichExclude),:);

numSubs = size(metadata,1);

% ------------------------------------------------------------------------------		
% Define covariates
% ------------------------------------------------------------------------------
Covs = {};
% Covs = {'Age','Gender'};
% Covs = {'Age','Gender','fdJenk_m'};
% Covs = {'Age','Gender','fdJenk_m','BDI','STAI_I'};

numCov = length(Covs);

% ------------------------------------------------------------------------------
% write out reduced metadata for subsequent analyses to read (e.g., spDCM)
% ------------------------------------------------------------------------------
extraStr = '';
baseoutdir = [projdir,'SecondLevel/SPM/Factorial/',WhichNoise,WhichSeed,extraStr,'/'];
% make baseoutdir
if exist(baseoutdir) == 0
	fprintf(1,'Initialising baseoutdir\n')
	mkdir(baseoutdir)
elseif exist(baseoutdir) == 7
	fprintf(1,'Cleaning and re-initialising baseoutdir\n')
	rmdir(baseoutdir,'s')
	mkdir(baseoutdir)
end
cd(baseoutdir)
save('metadata.mat','metadata','Covs')

% Also write out a .txt file listing subject IDs for MASSIVE_Job.sh scripts
fileID = fopen('ParticipantIDs.txt','w');
ParticipantIDs = metadata.ParticipantID;
fprintf(fileID,'%s\n',ParticipantIDs{:});
fclose(fileID);
clear ParticipantIDs

% ------------------------------------------------------------------------------
% SPM
% ------------------------------------------------------------------------------
for conNum = conNums
    outdir = [baseoutdir,'ROI_',num2str(conNum),'/'];
	mkdir(outdir)
	cd(outdir)

	%==========================================================================
	% Script
	%==========================================================================

	%-----------------------------------------------------------------------
	% Job configuration created by cfg_util (rev $Rev: 4252 $)
	%-----------------------------------------------------------------------

	fprintf(1,'Initialising SPM...\n')
	spm('defaults','fmri');
	spm_jobman('initcfg')
	spm_get_defaults('mask.thresh',-Inf)
	% spm_get_defaults('images.format','nii')

	% ------------------------------------------------------------------------------
	% Input con images to factor (hemi)
	% ------------------------------------------------------------------------------

	matlabbatch{1}.spm.stats.factorial_design.dir = {outdir};
	matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).name = 'Group';
	matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).levels = numGroups;
	matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).dept = 0;
	matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).variance = 1;
	matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).gmsca = 0;
	matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).ancova = 0;
	matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).name = 'Hemisphere';
	matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).levels = 2; % 2 levels corresponding to 1 = left and 2 = right
	matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).dept = 0;
	matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).variance = 1;
	matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).gmsca = 0;
	matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).ancova = 0;

	% ------------------------------------------------------------------------------
	% generate matrix that is then used to specify the cells and levels in the matlabbatch
	% ------------------------------------------------------------------------------
	WhichCell = 1:(2*numGroups);
	WhichCell = reshape(WhichCell,2,numGroups);
	WhichCell = WhichCell';

	% above matrix will be numGroups by Hemisphere (i.e., 3 groups by 2 hemispheres)
	% each element will correspond to the cell for the factorial design model
	% and the row and column indices of each element will correspond to the levels of the two factors
	% e.g.,
	% Three groups:
	% WhichCell =
	%      1     2
	%      3     4
	%      5     6
	% 
	% Two groups:
	% WhichCell =
	%      1     2
	%      3     4
	%
	% so WhichCell == 1 corresponds to element 1 which corresponds to row 1 and column 1
	% thus corresponding to cell 1 in the model which is group 1 hemisphere 1 (i.e., left)
	% WhichCell == 2 corresponds to element 2 which corresponds to row 1 column 2
	% thus corresponding to cell 2 in the model which is group 1 hemisphere 2 (i.e., right)
	% this will scale with the number of groups input to this script

	% Similarly, you can drop the right hemisphere easily. Like so:
	% WhichCell = 1:numGroups;
	% WhichCell = WhichCell';

	% WhichCell =
	%     1
	%     2
	% Dropping left hemisphere instead of right would requires changes to the below code.
	% Note, this is also what you'd use if you had already collapsed hemispheres during the first level model.

	numCells = numel(WhichCell);
	for i = 1:numCells
		% ------------------------------------------------------------------------------
		% find levels of ith cell
		% ------------------------------------------------------------------------------
		[group,hemi] = find(WhichCell == i); 
		matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(i).levels = [group,hemi];

		% ------------------------------------------------------------------------------
		% Specify which first level folder to import (i.e., left or right hemi)
		% ------------------------------------------------------------------------------
		if hemi == 1
			firstdir = [preprostr,WhichNoise,'FirstLevel/',WhichSeed,'_L/'];
			% firstdir = [preprostr,WhichNoise,'FirstLevel/',WhichSeed,'/'];
		elseif hemi == 2
			firstdir = [preprostr,WhichNoise,'FirstLevel/',WhichSeed,'_R/'];
		end
		% firstdir = [preprostr,WhichNoise,'FirstLevel/',WhichSeed,'_B/'];
		
		% ------------------------------------------------------------------------------
		% loop over all the subjects from group X
		% ------------------------------------------------------------------------------
		for j = 1:sum(metadata.Diagnosis == group)
			% input images for subject
			matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(i).scans{j,1} = [datadir,'/',metadata(metadata.Diagnosis == group,:).ParticipantID{j},firstdir,'/con_000',num2str(conNum),'.nii,1'];
			% matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(i).scans{j,1} = [datadir,'/',metadata(metadata.Diagnosis == group,:).ParticipantID{j},firstdir,'/con_000',num2str(conNum),'.img,1'];
		end
	end

	% ------------------------------------------------------------------------------
	% Covariates
	% ------------------------------------------------------------------------------
	if numCov > 0;
		% mean center age
		% metadata.Age = metadata.Age - mean(metadata.Age);

		for i = 1:numCov
			% initialise covariate
			c = [];

			% Get the covariate for all participants from table
			x = table2array(metadata(:,Covs{i}));
			for j = 1:numGroups
				if size(WhichCell,2) == 1
					y = x(metadata.Diagnosis == j);
				elseif size(WhichCell,2) == 2
					% for each group, duplicate to account for hemispheres being modeled separately in GLM
					y = [x(metadata.Diagnosis == j); x(metadata.Diagnosis == j)];
				end

				% store in c
				c = [c; y];
			end

			% input covariate for subject
			matlabbatch{1}.spm.stats.factorial_design.cov(i).c = c;
			matlabbatch{1}.spm.stats.factorial_design.cov(i).cname = Covs{i};
			matlabbatch{1}.spm.stats.factorial_design.cov(i).iCFI = 1;
			matlabbatch{1}.spm.stats.factorial_design.cov(i).iCC = 1;
		end
	else
		matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
	end

	matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
	matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
	matlabbatch{1}.spm.stats.factorial_design.masking.em = {''}; % '/media/lindenmp/ResProjects1/OCDPG/rfMRI_SecondLevel/mask.nii,1'
	matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
	matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
	matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

	% ------------------------------------------------------------------------------
	% Estimate
	% ------------------------------------------------------------------------------
	matlabbatch{2}.spm.stats.fmri_est.spmmat = {[outdir,'/SPM.mat']};
	matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

	% ------------------------------------------------------------------------------
	% Contrasts
	% ------------------------------------------------------------------------------
	matlabbatch{3}.spm.stats.con.spmmat = {[outdir,'/SPM.mat']};
	matlabbatch{3}.spm.stats.con.delete = 1;

	if size(WhichCell,2) == 1
		TheOnes = 1; 
		TheNegOnes = -1;
		TheZeros = 0;

		TheLOnes = 1; 
		TheROnes = 1; 
	elseif size(WhichCell,2) == 2
		TheOnes = [1 1]; 
		TheNegOnes = [-1 -1];
		TheZeros = [0 0];

		TheLOnes = [1 0]; 
		TheROnes = [0 1]; 
	end

	% Define contrasts
	if numGroups == 1
		% ------------------------------------------------------------------------------
		% Main effect of seed: ALL
		matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'Effect of Seed: All';
		matlabbatch{3}.spm.stats.con.consess{1}.tcon.convec = [TheOnes];
		matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';

		matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'Effect of Seed: All Left';
		matlabbatch{3}.spm.stats.con.consess{2}.tcon.convec = [TheLOnes];
		matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';

		matlabbatch{3}.spm.stats.con.consess{3}.tcon.name = 'Effect of Seed: All Right';
		matlabbatch{3}.spm.stats.con.consess{3}.tcon.convec = [TheROnes];
		matlabbatch{3}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
	elseif numGroups == 2
		% ------------------------------------------------------------------------------
		% Main effect of seed: ALL
		matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'Effect of Seed: All';
		matlabbatch{3}.spm.stats.con.consess{1}.tcon.convec = [TheOnes TheOnes];
		matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';

		matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'Effect of Seed: All Left';
		matlabbatch{3}.spm.stats.con.consess{2}.tcon.convec = [TheLOnes TheLOnes];
		matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';

		matlabbatch{3}.spm.stats.con.consess{3}.tcon.name = 'Effect of Seed: All Right';
		matlabbatch{3}.spm.stats.con.consess{3}.tcon.convec = [TheROnes TheROnes];
		matlabbatch{3}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
	elseif numGroups == 3
		% numContrasts = 2;
		% ------------------------------------------------------------------------------
		% Main effect of seed: ALL
		matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'Effect of Seed: All';
		matlabbatch{3}.spm.stats.con.consess{1}.tcon.convec = [TheOnes TheOnes TheOnes];
		matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';

		matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'Effect of Seed: All Left';
		matlabbatch{3}.spm.stats.con.consess{2}.tcon.convec = [TheLOnes TheLOnes TheLOnes];
		matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';

		matlabbatch{3}.spm.stats.con.consess{3}.tcon.name = 'Effect of Seed: All Right';
		matlabbatch{3}.spm.stats.con.consess{3}.tcon.convec = [TheROnes TheROnes TheROnes];
		matlabbatch{3}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
	end

	% ------------------------------------------------------------------------------
	% Run
	% ------------------------------------------------------------------------------
	spm_jobman('run',matlabbatch);
	clear matlabbatch
end
