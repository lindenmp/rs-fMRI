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
    	conNums = [2,3];
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
% Non-imaging data
% ------------------------------------------------------------------------------
data = readtable(sublist);

numGroups = length(unique(data.diagnosis));

% Generate parIDs var by group
parIDs = cell(1,numGroups);
for i = 1:numGroups
	parIDs{i} = data.participant_id(data.diagnosis == i);
end

numSubs = length(parIDs{2});

% ------------------------------------------------------------------------------
% Regressors
% ------------------------------------------------------------------------------
% predictor = 'StopSignal_SSRT'; Nothing at p<.001
predictor = 'DelayDiscounting_k';

% ------------------------------------------------------------------------------
% SPM
% ------------------------------------------------------------------------------
for conNum = conNums
    outdir = [projdir,'SecondLevelSPM/Regression/',WhichSeed,'/ROI_',num2str(conNum),'/'];

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

    %-----------------------------------------------------------------------
    % Job configuration created by cfg_util (rev $Rev: 4252 $)
    %-----------------------------------------------------------------------

    fprintf(1,'Initialising SPM...\n')
    spm('defaults','fmri');
    spm_jobman('initcfg')

    matlabbatch{1}.spm.stats.factorial_design.dir = {outdir};

    %-----------------------------------------------------------------------
    % Import scans
    %-----------------------------------------------------------------------

    for i = 1:numSubs
        % matlabbatch{1}.spm.stats.factorial_design.des.mreg.scans{i} = [datadir,'/',ParticipantIDs{i},firstdir,'con_000',num2str(conNum),'.img,1'];
        matlabbatch{1}.spm.stats.factorial_design.des.mreg.scans{i} = [datadir,'/',parIDs{2}{i},firstdir,'con_000',num2str(conNum),'.img,1'];
    end

    %-----------------------------------------------------------------------
    % Import predictors
    %-----------------------------------------------------------------------
    matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(1).c = data{data.diagnosis == 2,predictor};
    matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(1).cname = predictor;
    matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(1).iCC = 1; % 1 = overall mean centering, 5 = no mean centering

    % ------------------------------------------------------------------------------
    % Import covariates
    % ------------------------------------------------------------------------------    

    % for i = 1:numCov
    %     idx = strmatch(Covs{i},Bh,'exact');
        
    %     for j = 1:numSubs
    %         matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(i+numPred).c(j) = [B(j,idx)];
    %     end

    %     matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(i+numPred).cname = Covs{i};
    %     matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(i+numPred).iCC = 1;

    % end

    % ------------------------------------------------------------------------------
    % Design defaults
    % ------------------------------------------------------------------------------

    matlabbatch{1}.spm.stats.factorial_design.des.mreg.incint = 1;
    matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
    matlabbatch{1}.spm.stats.factorial_design.masking.im = 1; % 1 = implicit mask yes, 0 = implicit mask no
    matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
    matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

    % ------------------------------------------------------------------------------
    % Estimate
    % ------------------------------------------------------------------------------
    matlabbatch{2}.spm.stats.fmri_est.spmmat = {[outdir,'SPM.mat']};
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

    % ------------------------------------------------------------------------------
    % Contrasts
    % ------------------------------------------------------------------------------
    matlabbatch{3}.spm.stats.con.spmmat = {[outdir,'SPM.mat']};
    matlabbatch{3}.spm.stats.con.delete = 1;

    matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'Predictor 1+';
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.convec = [0 1];
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';

    matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'Predictor 1-';
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.convec = [0 -1];
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';

    % ------------------------------------------------------------------------------
    % Run
    % ------------------------------------------------------------------------------
    spm_jobman('run',matlabbatch);
    clear matlabbatch
end


