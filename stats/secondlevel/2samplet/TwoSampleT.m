% ------------------------------------------------------------------------------
% TwoSampleT.m
% ------------------------------------------------------------------------------
% This script will run a simple two sample t-test for the OCDPG dataset.
% This script was just used to check whether the ventral striatum - OFC effect
% commonly reported in the OCD literature was also present in BMH's OCD cohort.
% 
% Note, this script does not model hemisphere effects nor does it include any covariates
% ------------------------------------------------------------------------------
% Linden Parkes, Brain & Mental Health Laboratory, 2016
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

% Select which group to remove
% 2 = OCDs
% 3 = PGs
WhichGroup = 3;

% ------------------------------------------------------------------------------
% Set parcellation
% Note, code is not setup to process multiple parcellations concurrently.
% ------------------------------------------------------------------------------
switch WhichSeed
    case 'TriStri'
		conNum = 3;
    case 'DiMartino'
end

% ------------------------------------------------------------------------------
% Set project variables
% ------------------------------------------------------------------------------
switch WhichProject
    case 'OCDPG'
        projdir = '/gpfs/M2Home/lindenmp/Monash076/Linden/OCDPG/';
        sublist = [projdir,'OCDPGe.csv'];
        datadir = [projdir,'data/'];
        preprostr = '/rfMRI/prepro/';
		
		firstdir = [preprostr,'/',removeNoise,'FirstLevel_L_',WhichSeed,'/'];
		outdir = [projdir,'SecondLevelSPM/TwoSampleT/',removeNoise,'/',WhichSeed];
end

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

% ------------------------------------------------------------------------------
% Subject list
% ------------------------------------------------------------------------------
fileID = fopen(sublist);
metadata = textscan(fileID, '%s %u %u %s %s','HeaderLines',1, 'delimiter',',');
metadata{2} = double(metadata{2}); metadata{3} = double(metadata{3});
ParticipantIDs = metadata{1};
Group = metadata{2};

% Remove one group
ParticipantIDs(Group == WhichGroup) = [];
Group(Group == WhichGroup) = [];

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
