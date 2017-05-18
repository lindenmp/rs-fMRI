% ------------------------------------------------------------------------------
% TwoSampleT.m
% ------------------------------------------------------------------------------
% This script will run a simple two sample t-test for the OCDPG dataset.
% 
% Note, this script does not model hemisphere effects nor does it include any covariates
% ------------------------------------------------------------------------------
% Linden Parkes, Brain & Mental Health Laboratory, 2017
% ------------------------------------------------------------------------------
clear all; close all; clc

% ------------------------------------------------------------------------------
% Paths
% ------------------------------------------------------------------------------
% where the prepro scripts are
funcdir = '/gpfs/M2Home/projects/Monash076/Linden/scripts/rs-fMRI/func/';
addpath(funcdir)

% ------------------------------------------------------------------------------
% Set options
% ------------------------------------------------------------------------------
Projects = {'OCDPG'};
WhichProject = Projects{1};

removeNoise = 'sICA-AROMA+2P';

WhichSeed = 'TriStri';

% ------------------------------------------------------------------------------
% Set parcellation
% Note, code is not setup to process multiple parcellations concurrently.
% ------------------------------------------------------------------------------
switch WhichSeed
    case 'TriStri'
    	conNums = [1,2,3];
		% conNum = 3;
    case 'DiMartino'
end

% ------------------------------------------------------------------------------
% Set project variables
% ------------------------------------------------------------------------------
switch WhichProject
    case 'OCDPG'
        projdir = '/gpfs/M2Home/lindenmp/Monash076/Linden/OCDPG/';
        sublist = [projdir,'HCOCD.csv'];
        datadir = [projdir,'data/'];
        preprostr = '/rfMRI/prepro/';
		
		firstdir = [preprostr,'/',removeNoise,'/FirstLevel_L_',WhichSeed,'/'];
end

% ------------------------------------------------------------------------------
% Subject list
% ------------------------------------------------------------------------------
fileID = fopen(sublist);
metadata = textscan(fileID, ['%s %f %f %s %s' repmat(' %f', 1, 41)],'HeaderLines',1, 'delimiter',',');
ParticipantIDs = metadata{1};
Group = metadata{2};

numSubs = length(ParticipantIDs);
numGroups = numel(unique(Group));
    
% Generate parIDs var by group
parIDs = cell(1,numGroups);
for i = 1:numGroups
	parIDs{i} = ParticipantIDs(Group == i);
end

% ------------------------------------------------------------------------------
% SPM
% ------------------------------------------------------------------------------
for conNum = conNums
	outdir = [projdir,'SecondLevelSPM/TwoSampleT/',removeNoise,'/',WhichSeed,'/contrast_',num2str(conNum),'/'];

	% make outdir
	if exist(outdir) == 0
		fprintf(1,'Initialising outdir\n')
		mkdir(outdir)
	elseif exist(outdir) == 7
		fprintf(1,'Cleaning and re-initialising outdir\n')
		rmdir(outdir,'s')
		mkdir(outdir)
	end

	cd(outdir)

	fprintf(1,'Initialising SPM...\n')
	spm('defaults','fmri');
	spm_jobman('initcfg')

	matlabbatch{1}.spm.stats.factorial_design.dir = {outdir};

	for i = 1:length(parIDs{1})
		matlabbatch{1}.spm.stats.factorial_design.des.t2.scans1{i,1} = [datadir,'/',parIDs{1}{i},firstdir,'con_000',num2str(conNum),'.img,1'];
	end

	for i = 1:length(parIDs{2})
		matlabbatch{1}.spm.stats.factorial_design.des.t2.scans2{i,1} = [datadir,'/',parIDs{2}{i},firstdir,'con_000',num2str(conNum),'.img,1'];
	end

	matlabbatch{1}.spm.stats.factorial_design.des.t2.dept = 0;
	matlabbatch{1}.spm.stats.factorial_design.des.t2.variance = 1;
	matlabbatch{1}.spm.stats.factorial_design.des.t2.gmsca = 0;
	matlabbatch{1}.spm.stats.factorial_design.des.t2.ancova = 0;
	matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
	matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
	matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
	matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
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

		% ------------------------------------------------------------------------------
		% T contrasts
		matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'HC > Patients';
		matlabbatch{3}.spm.stats.con.consess{1}.tcon.convec = [1 -1];
		matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';

		matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'Patients > HC';
		matlabbatch{3}.spm.stats.con.consess{2}.tcon.convec = [-1 1];
		matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';

		% ------------------------------------------------------------------------------
		% T contrasts
		matlabbatch{3}.spm.stats.con.consess{3}.tcon.name = 'HC main effect';
		matlabbatch{3}.spm.stats.con.consess{3}.tcon.convec = [1 0];
		matlabbatch{3}.spm.stats.con.consess{3}.tcon.sessrep = 'none';

		matlabbatch{3}.spm.stats.con.consess{4}.tcon.name = 'Patients main effect';
		matlabbatch{3}.spm.stats.con.consess{4}.tcon.convec = [0 1];
		matlabbatch{3}.spm.stats.con.consess{4}.tcon.sessrep = 'none';

	% ------------------------------------------------------------------------------
	% Run
	% ------------------------------------------------------------------------------
	spm_jobman('run',matlabbatch);
	clear matlabbatch
end
