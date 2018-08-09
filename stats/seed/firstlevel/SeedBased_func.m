%% SeedBased: 
function [] = SeedBased_func(WhichProject,ParticipantID,removeNoise,WhichSeed)
    % ------------------------------------------------------------------------------
    % SeedBased.m
    % ------------------------------------------------------------------------------
    % This script will run the first level analysis.
    % 
    % ------------------------------------------------------------------------------
    % Linden Parkes, Brain & Mental Health Laboratory, 2016
    % ------------------------------------------------------------------------------

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

    fprintf(1, 'Running first level analysis using %s ROIs\n',WhichSeed);

    % ------------------------------------------------------------------------------
    % Set parcellation
    % Note, code is not setup to process multiple parcellations concurrently.
    % ------------------------------------------------------------------------------
    switch WhichSeed
        case 'TriStri'
            Seed = 3;
            idx = [1,3; 4,6; 1,6]; % left / right / bilat
        case 'DiMartino'
            Seed = 4;
            idx = [1,6; 7,12; 1,12];
    end

    % ------------------------------------------------------------------------------
    % Set project variables
    % ------------------------------------------------------------------------------
    switch WhichProject
        case 'OCDPG_DCM'
            datadir = [parentdir_scratch,'ResProjects/rfMRI_DCM/OCDPG/data/'];
            preprostr = '/func/prepro/';

            TR = 2.5;
            N = 189;
            N = N - 4;
            numSlices = 44;
            refSlice = numSlices - 1; % use for interleaved order

            % Name of preprocessed functional data file
            EPI = 'epi_prepro.nii';
            brainMask = 'brain_mask.nii';
        case 'GenCog'
            datadir = [parentdir_scratch,'ResProjects/rfMRI_DCM/GenCog/data/'];
            preprostr = '/func/prepro/';        

            TR = 0.754;
            N = 616;
            numSlices = 16;
            refSlice = 8;

            % Name of preprocessed functional data file
            EPI = 'epi_prepro.nii';
            brainMask = 'MNI152_T1_2mm_brain_mask.nii';
        case 'UCLA'
            datadir = [parentdir_scratch,'ResProjects/rfMRI_denoise/UCLA/data/'];
            preprostr = '/func/prepro/';        

            TR = 2;
            N = 152;
            N = N - 4;
            numSlices = 34;
            refSlice = numSlices - 1; % use for interleaved order

            % Name of preprocessed functional data file
            EPI = 'epi_prepro.nii';
            brainMask = 'brain_mask.nii';
    end

    % ------------------------------------------------------------------------------
    % Script
    % ------------------------------------------------------------------------------
    fprintf(1,'Processing subject %s\n',ParticipantID)

    preprodir = [datadir,ParticipantID,preprostr];
    cleandir = [preprodir,removeNoise,'/'];

    cd(cleandir)
    if exist(EPI) == 0
        if exist([EPI,'.gz']) == 2
            runDecompEPI = 1;
        elseif exist([EPI,'.gz']) == 0
            fprintf(1, 'Warning: EPI not found!\n');
        end
    elseif exist(EPI) == 2;
        runDecompEPI = 0;
    end

    if runDecompEPI == 1
        fprintf(1, '\t\t Decompressing epi...\n');
        gunzip([EPI,'.gz'])
        delete([EPI,'.gz'])
    end

    cd(preprodir)
    if exist(brainMask) == 0
        if exist([brainMask,'.gz']) == 2
            runDecompbrainMask = 1;
        elseif exist([brainMask,'.gz']) == 0
            fprintf(1, 'Warning: brainMask not found!\n');
        end
    elseif exist(brainMask) == 2;
        runDecompbrainMask = 0;
    end

    if runDecompbrainMask == 1
        fprintf(1, '\t\t Decompressing brainMask...\n');
        gunzip([brainMask,'.gz'])
        delete([brainMask,'.gz'])
    end

    % ------------------------------------------------------------------------------
    % Get time series
    % ------------------------------------------------------------------------------
    if ismember('TriStri',WhichSeed,'rows') | ismember('DiMartino',WhichSeed,'rows')
        cd(cleandir)
        load('cfg.mat');
        roiTS = cfg.roiTS{Seed};
    elseif ismember('GSR',WhichSeed,'rows')
        roiTS = dlmread([preprodir,'6P+2P+GSR/gsTS.txt']);
    end

    % ------------------------------------------------------------------------------
    % Estimate first level
    % ------------------------------------------------------------------------------
    workdir = [cleandir,'FirstLevel/'];
    if ismember('TriStri',WhichSeed,'rows') | ismember('DiMartino',WhichSeed,'rows')
        firstDirs = {[workdir,WhichSeed,'_L/'],[workdir,WhichSeed,'_R/'],[workdir,WhichSeed,'_B/']};
    elseif ismember('GSR',WhichSeed,'rows')
        firstDirs = {[workdir,WhichSeed,'/']};
    end

    for j = 1:length(firstDirs)
        firstdir = firstDirs{j};

        % Setup output directory
        if exist(firstdir) == 0
            fprintf(1,'\t\tInitialising firstdir\n')
            mkdir(firstdir)
        elseif exist(firstdir) == 7
            fprintf(1,'\t\tCleaning and re-initialising firstdir\n')
            rmdir(firstdir,'s')
            mkdir(firstdir)
        end

        if ismember('TriStri',WhichSeed,'rows') | ismember('DiMartino',WhichSeed,'rows')
            R = roiTS(:,idx(j,1):idx(j,2));
        elseif ismember('GSR',WhichSeed,'rows')
            R = roiTS;
        end
        
        save spm_regs R
        clear R
        movefile('spm_regs.mat',firstdir)

        FirstLevelGLM(firstdir,'scans',TR,numSlices,refSlice,[cleandir,EPI],N,[firstdir,'spm_regs.mat'],[preprodir,brainMask]);
        
        if j == 3
            FirstLevelContrasts([WhichSeed,'_B'],[firstdir,'SPM.mat']);
        else
            FirstLevelContrasts(WhichSeed,[firstdir,'SPM.mat']);
        end
    end
    fprintf(1, 'First level analysis done \n');
end
