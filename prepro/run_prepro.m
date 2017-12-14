%% run_prepro:
% Copyright (C) 2017, Linden Parkes <lindenparkes@gmail.com>,

function [] = run_prepro(WhichProject,WhichSessScan,subject,smoothing,discard,slicetime,despike,detr,intnorm,bandpass,demean)
    cfg.WhichProject = WhichProject;
    cfg.WhichSessScan = WhichSessScan;
    cfg.subject = subject;
    % preprocessing options
    if nargin < 5
        cfg.discard = 1;
    else
        cfg.discard = discard;
    end

    if nargin < 6
        cfg.slicetime = 1;
    else
        cfg.slicetime = slicetime;
    end

    if nargin < 7
        cfg.despike = 0;
    else
        cfg.despike = despike;
    end

    if nargin < 8
        cfg.detr = 1;
    else
        cfg.detr = detr;
    end

    if nargin < 9
        cfg.intnorm = 1;
    else
        cfg.intnorm = intnorm;
    end

    if nargin < 10
        cfg.runBandpass = 1;
    else
        cfg.runBandpass = bandpass;
    end

    if nargin < 11
        cfg.demean = 1;
    else
        cfg.demean = demean;
    end

    % ------------------------------------------------------------------------------
    % Store date and time
    % ------------------------------------------------------------------------------
    cfg.DateTime = datetime('now');

    % ------------------------------------------------------------------------------
    % Parent dir
    % ------------------------------------------------------------------------------
    cfg.parentdir = '/home/lindenmp/kg98/Linden/';
    cfg.parentdir_scratch = '/home/lindenmp/kg98_scratch/Linden/';

    % ------------------------------------------------------------------------------
    % Add paths - edit this section
    % ------------------------------------------------------------------------------
    % where the prepro scripts are
    cfg.scriptdir = [cfg.parentdir,'Scripts/rs-fMRI/prepro/'];
    addpath(cfg.scriptdir)
    cfg.funcdir = [cfg.parentdir,'Scripts/rs-fMRI/func/'];
    addpath(cfg.funcdir)

    % where spm is
    cfg.spmdir = '/usr/local/spm8/matlab2015b.r6685/';
    addpath(cfg.spmdir)

    % set FSL environments
    cfg.fsldir = '/usr/local/fsl/5.0.9/fsl/bin/';
    setenv('FSLDIR',cfg.fsldir(1:end-4));
    setenv('FSLOUTPUTTYPE','NIFTI');
    setenv('LD_LIBRARY_PATH',[getenv('PATH'),getenv('LD_LIBRARY_PATH'),':/usr/lib/fsl/5.0'])

    % where ICA-AROMA scripts are
    cfg.scriptdir_ICA = [cfg.parentdir,'Scripts/Tools/ICA-AROMA/'];

    % ANTs
    cfg.antsdir = '/usr/local/ants/1.9.v4/bin/';
    setenv('ANTSPATH',cfg.antsdir);

    % Directory to AFNI functions
    cfg.afnidir = '/usr/local/afni/16.2.16/';

    % ------------------------------------------------------------------------------
    % Set project settings and parameters
    % Use WhichProject if you're juggling multiple datasets
            % cfg.order = [1:1:cfg.numSlices]; % ascending
            % cfg.order = [cfg.numSlices:-1:1]; % descending
            % cfg.order = [1:2:cfg.numSlices,2:2:cfg.numSlices]; % interleaved
            % cfg.order = [2:2:cfg.numSlices, 1:2:cfg.numSlices-1]; %interleave alt
            % Reference slice for slice timing acquisition. See SlicetimeEPI.m for guidance on how to define.
            % cfg.refSlice = round(cfg.numSlices/2); % use for sequential order (e.g., ascending or descending)
            % cfg.refSlice = cfg.numSlices-1; % use for interleaved order
            % cfg.refSlice = cfg.numSlices; % use for interleaved alt order
    % ------------------------------------------------------------------------------
    switch cfg.WhichProject
        case 'OCDPG'
            % Where the subjects' directories are
            cfg.datadir = [cfg.parentdir_scratch,'ResProjects/rfMRI_denoise/OCDPG/data/'];

            switch cfg.WhichSessScan
                case 'Sess1_Scan1'
                    % where the unprocessed EPI 4d file is
                    cfg.rawdir = [cfg.datadir,cfg.subject,'/func/'];
                    % Directory where the t1 is
                    cfg.t1dir = [cfg.datadir,cfg.subject,'/anat/'];
            end

            % where the processed epi 4d files will be output to from prepro_base
            cfg.preprodir = [cfg.rawdir,'prepro/'];

            % file name of EPI 4d file
            cfg.EPI = [cfg.subject,'_task-rest_bold.nii'];
            % name of t1 file.
            cfg.t1name = [cfg.subject,'_T1w.nii'];

            % the path and filename of the template in MNI space to which everything
            % will be normalized
            cfg.mni_template = [cfg.spmdir,'templates/T1.nii'];

            % preprocessing settings
            % length of time series (no. vols)
            cfg.N = 189;
            % Repetition time of acquistion in secs
            cfg.TR = 2.5;
            % Number of slices in EPI volumes.
            cfg.numSlices = 44;
            % Vector defining acquisition order of EPI slices (necessary for slice-timing correction.
            % See help of slicetime_epis.m for guidance on how to define)
            cfg.order = [1:2:cfg.numSlices,2:2:cfg.numSlices]; % interleaved
            % Reference slice for slice timing acquisition. See SlicetimeEPI.m for guidance on how to define.
            cfg.refSlice = cfg.numSlices-1; % use for interleaved order

            % Scalar value indicating spatial smoothing kernal size in mm.
            cfg.kernel = 6;
            % Low-pass cut-off for bandpass filter in Hz (e.g., .08)
            cfg.LowPass = 0.08;
            % Hi-pass cut-off for bandpass filter in Hz (e.g., .008)
            cfg.HighPass = 0.008;
        case 'UCLA'
            % Where the subjects' directories are
            cfg.datadir = [cfg.parentdir_scratch,'ResProjects/rfMRI_denoise/UCLA/data/'];

            switch cfg.WhichSessScan
                case 'Sess1_Scan1'
                    % where the unprocessed EPI 4d file is
                    cfg.rawdir = [cfg.datadir,cfg.subject,'/func/'];
                    % Directory where the t1 is
                    cfg.t1dir = [cfg.datadir,cfg.subject,'/anat/'];
            end

            % where the processed epi 4d files will be output to from prepro_base
            cfg.preprodir = [cfg.rawdir,'prepro/'];

            % file name of EPI 4d file
            cfg.EPI = [cfg.subject,'_task-rest_bold.nii'];
            % name of t1 file.
            cfg.t1name = [cfg.subject,'_T1w.nii'];

            % the path and filename of the template in MNI space to which everything
            % will be normalized
            cfg.mni_template = [cfg.spmdir,'templates/T1.nii'];

            % preprocessing settings
            % length of time series (no. vols)
            cfg.N = 152;
            % Repetition time of acquistion in secs
            cfg.TR = 2;
            % Number of slices in EPI volumes.
            cfg.numSlices = 34;
            % Vector defining acquisition order of EPI slices (necessary for slice-timing correction.
            % See help of slicetime_epis.m for guidance on how to define)
            cfg.order = [1:2:cfg.numSlices,2:2:cfg.numSlices]; % interleaved
            % Reference slice for slice timing acquisition. See SlicetimeEPI.m for guidance on how to define.
            cfg.refSlice = cfg.numSlices-1; % use for interleaved order

            % Scalar value indicating spatial smoothing kernal size in mm.
            cfg.kernel = 6;
            % Low-pass cut-off for bandpass filter in Hz (e.g., .08)
            cfg.LowPass = 0.08;
            % Hi-pass cut-off for bandpass filter in Hz (e.g., .008)
            cfg.HighPass = 0.008;
        case 'NYU_2'
            % Where the subjects' directories are
            cfg.datadir = [cfg.parentdir_scratch,'ResProjects/rfMRI_denoise/NYU_2/data/'];

            % where the unprocessed EPI 4d file is
            % cfg.WhichSessScan = 'Sess1_Scan1';
            switch cfg.WhichSessScan
                case 'Sess1_Scan1'
                    % where the unprocessed EPI 4d file is
                    cfg.rawdir = [cfg.datadir,cfg.subject,'/session_1/func_1/'];
                    % Directory where the t1 is
                    cfg.t1dir = [cfg.datadir,cfg.subject,'/session_1/anat_1/'];
                case 'Sess1_Scan2'
                    cfg.rawdir = [cfg.datadir,cfg.subject,'/session_1/func_2/'];
                    cfg.t1dir = [cfg.datadir,cfg.subject,'/session_1/anat_1/'];
                case 'Sess2_Scan1'
                    cfg.rawdir = [cfg.datadir,cfg.subject,'/session_2/func_1/'];
                    cfg.t1dir = [cfg.datadir,cfg.subject,'/session_2/anat_1/'];
            end

            % where the processed epi 4d files will be output to from prepro_base
            cfg.preprodir = [cfg.rawdir,'prepro/'];

            % file name of EPI 4d file
            cfg.EPI = [cfg.subject,'_task-rest_bold.nii'];
            % name of t1 file.
            cfg.t1name = [cfg.subject,'_T1w.nii'];

            % the path and filename of the template in MNI space to which everything
            % will be normalized
            cfg.mni_template = [cfg.spmdir,'templates/T1.nii'];

            % preprocessing settings
            % length of time series (no. vols)
            cfg.N = 180;
            % Repetition time of acquistion in secs
            cfg.TR = 2;
            % Number of slices in EPI volumes.
            cfg.numSlices = 33;
            % Vector defining acquisition order of EPI slices (necessary for slice-timing correction.
            % See help of slicetime_epis.m for guidance on how to define)
            cfg.order = [1:2:cfg.numSlices,2:2:cfg.numSlices]; % interleaved
            % Reference slice for slice timing acquisition. See SlicetimeEPI.m for guidance on how to define.
            cfg.refSlice = cfg.numSlices-1; % use for interleaved order

            % Scalar value indicating spatial smoothing kernal size in mm.
            cfg.kernel = 6;
            % Low-pass cut-off for bandpass filter in Hz (e.g., .08)
            cfg.LowPass = 0.08;
            % Hi-pass cut-off for bandpass filter in Hz (e.g., .008)
            cfg.HighPass = 0.008;
        case 'COBRE'
            % Where the subjects' directories are
            cfg.datadir = [cfg.parentdir_scratch,'ResProjects/rfMRI_denoise/COBRE/data/'];

            switch cfg.WhichSessScan
                case 'Sess1_Scan1'
                    % where the unprocessed EPI 4d file is
                    cfg.rawdir = [cfg.datadir,cfg.subject,'/func/'];
                    % Directory where the t1 is
                    cfg.t1dir = [cfg.datadir,cfg.subject,'/anat/'];
            end

            % where the processed epi 4d files will be output to from prepro_base
            cfg.preprodir = [cfg.rawdir,'prepro/'];

            % file name of EPI 4d file
            cfg.EPI = [cfg.subject,'_task-rest_bold.nii'];
            % name of t1 file.
            cfg.t1name = [cfg.subject,'_T1w.nii'];

            % the path and filename of the template in MNI space to which everything
            % will be normalized
            cfg.mni_template = [cfg.spmdir,'templates/T1.nii'];

            % preprocessing settings
            % length of time series (no. vols)
            cfg.N = 150;
            % Repetition time of acquistion in secs
            cfg.TR = 2;
            % Number of slices in EPI volumes.
            cfg.numSlices = 33;
            % Vector defining acquisition order of EPI slices (necessary for slice-timing correction.
            % See help of slicetime_epis.m for guidance on how to define)
            cfg.order = [1:2:cfg.numSlices,2:2:cfg.numSlices]; % interleaved
            % Reference slice for slice timing acquisition. See SlicetimeEPI.m for guidance on how to define.
            cfg.refSlice = cfg.numSlices-1; % use for interleaved order

            % Scalar value indicating spatial smoothing kernal size in mm.
            cfg.kernel = 6;
            % Low-pass cut-off for bandpass filter in Hz (e.g., .08)
            cfg.LowPass = 0.08;
            % Hi-pass cut-off for bandpass filter in Hz (e.g., .008)
            cfg.HighPass = 0.008;
        case 'Beijing_Zang'
            % Where the subjects' directories are
            cfg.datadir = [cfg.parentdir_scratch,'ResProjects/rfMRI_denoise/Beijing_Zang/data/'];

            switch cfg.WhichSessScan
                case 'Sess1_Scan1'
                    % where the unprocessed EPI 4d file is
                    cfg.rawdir = [cfg.datadir,cfg.subject,'/func/'];
                    % Directory where the t1 is
                    cfg.t1dir = [cfg.datadir,cfg.subject,'/anat/'];
            end

            % where the processed epi 4d files will be output to from prepro_base
            cfg.preprodir = [cfg.rawdir,'prepro/'];

            % file name of EPI 4d file
            cfg.EPI = [cfg.subject,'_task-rest_bold.nii'];
            % name of t1 file.
            cfg.t1name = [cfg.subject,'_T1w.nii'];

            % the path and filename of the template in MNI space to which everything
            % will be normalized
            cfg.mni_template = [cfg.spmdir,'templates/T1.nii'];

            % preprocessing settings
            % length of time series (no. vols)
            cfg.N = 225;
            % Repetition time of acquistion in secs
            cfg.TR = 2;
            % Number of slices in EPI volumes.
            cfg.numSlices = 33;
            % Vector defining acquisition order of EPI slices (necessary for slice-timing correction.
            % See help of slicetime_epis.m for guidance on how to define)
            cfg.order = [1:2:cfg.numSlices,2:2:cfg.numSlices]; % interleaved
            % Reference slice for slice timing acquisition. See SlicetimeEPI.m for guidance on how to define.
            cfg.refSlice = cfg.numSlices-1; % use for interleaved order

            % Scalar value indicating spatial smoothing kernal size in mm.
            cfg.kernel = 6;
            % Low-pass cut-off for bandpass filter in Hz (e.g., .08)
            cfg.LowPass = 0.08;
            % Hi-pass cut-off for bandpass filter in Hz (e.g., .008)
            cfg.HighPass = 0.008;
        case 'COBRE_HCTSA'
            % Where the subjects' directories are
            cfg.datadir = [cfg.parentdir,'ResProjects/SCZ_HCTSA/COBRE/data/'];

            switch cfg.WhichSessScan
                case 'Sess1_Scan1'
                    % where the unprocessed EPI 4d file is
                    cfg.rawdir = [cfg.datadir,cfg.subject,'/session_1/rest_1/'];
                    % Directory where the t1 is
                    cfg.t1dir = [cfg.datadir,cfg.subject,'/session_1/anat_1/'];
            end

            % where the processed epi 4d files will be output to from prepro_base
            cfg.preprodir = [cfg.rawdir,'prepro/'];

            % file name of EPI 4d file
            cfg.EPI = 'rest.nii';
            % name of t1 file.
            cfg.t1name = 'mprage.nii';

            % the path and filename of the template in MNI space to which everything
            % will be normalized
            cfg.mni_template = [cfg.spmdir,'templates/T1.nii'];

            % preprocessing settings
            % length of time series (no. vols)
            cfg.N = 150;
            % Repetition time of acquistion in secs
            cfg.TR = 2;
            % Number of slices in EPI volumes.
            cfg.numSlices = 33;
            % Vector defining acquisition order of EPI slices (necessary for slice-timing correction.
            % See help of slicetime_epis.m for guidance on how to define)
            cfg.order = [1:2:cfg.numSlices,2:2:cfg.numSlices]; % interleaved
            % Reference slice for slice timing acquisition. See SlicetimeEPI.m for guidance on how to define.
            cfg.refSlice = cfg.numSlices-1; % use for interleaved order

            % Scalar value indicating spatial smoothing kernal size in mm.
            cfg.kernel = 8;
            % Low-pass cut-off for bandpass filter in Hz (e.g., .08)
            cfg.LowPass = 10^6;
            % Hi-pass cut-off for bandpass filter in Hz (e.g., .008)
            cfg.HighPass = 0.008;
        case 'UCLA_HCTSA'
            % Where the subjects' directories are
            cfg.datadir = [cfg.parentdir,'ResProjects/SCZ_HCTSA/UCLA/data/'];

            switch cfg.WhichSessScan
                case 'Sess1_Scan1'
                    % where the unprocessed EPI 4d file is
                    cfg.rawdir = [cfg.datadir,cfg.subject,'/func/'];
                    % Directory where the t1 is
                    cfg.t1dir = [cfg.datadir,cfg.subject,'/anat/'];
            end

            % where the processed epi 4d files will be output to from prepro_base
            cfg.preprodir = [cfg.rawdir,'prepro/'];

            % file name of EPI 4d file
            cfg.EPI = [cfg.subject,'_task-rest_bold.nii'];
            % name of t1 file.
            cfg.t1name = [cfg.subject,'_T1w.nii'];

            % the path and filename of the template in MNI space to which everything
            % will be normalized
            cfg.mni_template = [cfg.spmdir,'templates/T1.nii'];

            % preprocessing settings
            % length of time series (no. vols)
            cfg.N = 152;
            % Repetition time of acquistion in secs
            cfg.TR = 2;
            % Number of slices in EPI volumes.
            cfg.numSlices = 34;
            % Vector defining acquisition order of EPI slices (necessary for slice-timing correction.
            % See help of slicetime_epis.m for guidance on how to define)
            cfg.order = [1:2:cfg.numSlices,2:2:cfg.numSlices]; % interleaved
            % Reference slice for slice timing acquisition. See SlicetimeEPI.m for guidance on how to define.
            cfg.refSlice = cfg.numSlices-1; % use for interleaved order

            % Scalar value indicating spatial smoothing kernal size in mm.
            cfg.kernel = 8;
            % Low-pass cut-off for bandpass filter in Hz (e.g., .08)
            cfg.LowPass = 10^6;
            % Hi-pass cut-off for bandpass filter in Hz (e.g., .008)
            cfg.HighPass = 0.008;
        case 'NAMIC_HCTSA'
            % Where the subjects' directories are
            cfg.datadir = [cfg.parentdir,'ResProjects/SCZ_HCTSA/NAMIC/data/'];

            switch cfg.WhichSessScan
                case 'Sess1_Scan1'
                    % where the unprocessed EPI 4d file is
                    cfg.rawdir = [cfg.datadir,cfg.subject,'/func/'];
                    % Directory where the t1 is
                    cfg.t1dir = [cfg.datadir,cfg.subject,'/anat/'];
            end

            % where the processed epi 4d files will be output to from prepro_base
            cfg.preprodir = [cfg.rawdir,'prepro/'];

            % file name of EPI 4d file
            cfg.EPI = [cfg.subject,'_task-rest_bold.nii'];
            % name of t1 file.
            cfg.t1name = [cfg.subject,'_T1w.nii'];

            % the path and filename of the template in MNI space to which everything
            % will be normalized
            cfg.mni_template = [cfg.spmdir,'templates/T1.nii'];

            % preprocessing settings
            % length of time series (no. vols)
            cfg.N = 200;
            % Repetition time of acquistion in secs
            cfg.TR = 3;
            % Number of slices in EPI volumes.
            cfg.numSlices = 39;
            % Vector defining acquisition order of EPI slices (necessary for slice-timing correction.
            % See help of slicetime_epis.m for guidance on how to define)
            cfg.order = [1:2:cfg.numSlices,2:2:cfg.numSlices]; % interleaved
            % Reference slice for slice timing acquisition. See SlicetimeEPI.m for guidance on how to define.
            cfg.refSlice = cfg.numSlices-1; % use for interleaved order

            % Scalar value indicating spatial smoothing kernal size in mm.
            cfg.kernel = 8;
            % Low-pass cut-off for bandpass filter in Hz (e.g., .08)
            cfg.LowPass = 10^6;
            % Hi-pass cut-off for bandpass filter in Hz (e.g., .008)
            cfg.HighPass = 0.008;
        case 'conc_TMS_FMRI'
            % Where the subjects' directories are
            cfg.datadir = '/home/lindenmp/kg98/Morrowj/PROJECTS/conc_TMS_FMRI/rawdata/';

            switch cfg.WhichSessScan
                case 'FUp_Scan1'
                    % where the unprocessed EPI 4d file is
                    cfg.rawdir = [cfg.datadir,cfg.subject,'/func/'];
                    cfg.preprodir = [cfg.rawdir,'prepro_FUp1/'];
                    cfg.EPI = [cfg.subject,'_task-FUP1_bold.nii'];

                    % Directory where the t1 is
                    cfg.t1dir = [cfg.datadir,cfg.subject,'/anat_FUp/']; % DLPFC, SFG and SHAM
                    cfg.t1name = [cfg.subject,'_FUP_T1w.nii'];
                case 'FUp_Scan2'
                    % where the unprocessed EPI 4d file is
                    cfg.rawdir = [cfg.datadir,cfg.subject,'/func/'];
                    cfg.preprodir = [cfg.rawdir,'prepro_FUp2/'];
                    cfg.EPI = [cfg.subject,'_task-FUP2_bold.nii'];

                    % Directory where the t1 is
                    cfg.t1dir = [cfg.datadir,cfg.subject,'/anat_FUp/']; % TPJ
                    cfg.t1name = [cfg.subject,'_FUP_T1w.nii'];
                case 'FDown_Scan1'
                    % where the unprocessed EPI 4d file is
                    cfg.rawdir = [cfg.datadir,cfg.subject,'/func/'];
                    cfg.preprodir = [cfg.rawdir,'prepro_FDown1/'];
                    cfg.EPI = [cfg.subject,'_task-FDOWN1_bold.nii'];

                    % Directory where the t1 is
                    cfg.t1dir = [cfg.datadir,cfg.subject,'/anat_FDown/']; % TPJ
                    cfg.t1name = [cfg.subject,'_FDOWN_T1w.nii'];
                    cfg.t14norm = [cfg.subject,'_FUP_T1w.nii'];
                case 'FDown_Scan2'
                    % where the unprocessed EPI 4d file is
                    cfg.rawdir = [cfg.datadir,cfg.subject,'/func/'];
                    cfg.preprodir = [cfg.rawdir,'prepro_FDown2/'];
                    cfg.EPI = [cfg.subject,'_task-FDOWN2_bold.nii'];

                    % Directory where the t1 is
                    cfg.t1dir = [cfg.datadir,cfg.subject,'/anat_FDown/']; % TPJ
                    cfg.t1name = [cfg.subject,'_FDOWN_T1w.nii'];
                    cfg.t14norm = [cfg.subject,'_FUP_T1w.nii'];
            end

            % the path and filename of the template in MNI space to which everything
            % will be normalized
            cfg.mni_template = [cfg.spmdir,'templates/T1.nii'];

            % preprocessing settings
            % length of time series (no. vols)
            cfg.N = 186;
            % Repetition time of acquistion in secs
            cfg.TR = 2.76;
            % Number of slices in EPI volumes.
            cfg.numSlices = 44;
            % Vector defining acquisition order of EPI slices (necessary for slice-timing correction.
            % See help of slicetime_epis.m for guidance on how to define)
            cfg.order = [1:2:cfg.numSlices,2:2:cfg.numSlices]; % interleaved
            % Reference slice for slice timing acquisition. See SlicetimeEPI.m for guidance on how to define.
            cfg.refSlice = cfg.numSlices-1; % use for interleaved order

            % Scalar value indicating spatial smoothing kernal size in mm.
            cfg.kernel = 6;
            % Low-pass cut-off for bandpass filter in Hz (e.g., .08)
            cfg.LowPass = 0.08;
            % Hi-pass cut-off for bandpass filter in Hz (e.g., .008)
            cfg.HighPass = 0.008;
        case 'Activebrains_pre'
            % Where the subjects' directories are
            cfg.datadir = '/scratch/kg98/IreneEsteban/Activebrains_pre/data/';

            switch cfg.WhichSessScan
                case 'Sess1_Scan1'
                    % where the unprocessed EPI 4d file is
                    cfg.rawdir = [cfg.datadir,cfg.subject,'/func/'];
                    % Directory where the t1 is
                    cfg.t1dir = [cfg.datadir,cfg.subject,'/anat/'];
            end

            % where the processed epi 4d files will be output to from prepro_base
            cfg.preprodir = [cfg.rawdir,'prepro/'];

            % file name of EPI 4d file
            cfg.EPI = [cfg.subject,'_task-rest_bold.nii'];
            % name of t1 file.
            cfg.t1name = [cfg.subject,'_T1w.nii'];

            % the path and filename of the template in MNI space to which everything
            % will be normalized
            cfg.mni_template = [cfg.spmdir,'templates/T1.nii'];

            % preprocessing settings
            % length of time series (no. vols)
            cfg.N = 160;
            % Repetition time of acquistion in secs
            cfg.TR = 2;
            % Number of slices in EPI volumes.
            cfg.numSlices = 35;
            % Vector defining acquisition order of EPI slices (necessary for slice-timing correction.
            % See help of slicetime_epis.m for guidance on how to define)
            cfg.order = [cfg.numSlices:-1:1]; % descending
            % Reference slice for slice timing acquisition. See SlicetimeEPI.m for guidance on how to define.
            cfg.refSlice = round(cfg.numSlices/2); % use for sequential order (e.g., ascending or descending)

            % Scalar value indicating spatial smoothing kernal size in mm.
            cfg.kernel = 6;
            % Low-pass cut-off for bandpass filter in Hz (e.g., .08)
            cfg.LowPass = 0.08;
            % Hi-pass cut-off for bandpass filter in Hz (e.g., .008)
            cfg.HighPass = 0.008;
        case 'TBS_fMRI'
            % Where the subjects' directories are
            cfg.datadir = [cfg.parentdir_scratch,'ResProjects/TBS_fMRI/data/'];

            switch cfg.WhichSessScan
                case 'Sess1_Scan1'
                    % where the unprocessed EPI 4d file is
                    cfg.rawdir = [cfg.datadir,cfg.subject,'/func_1/'];
                    % Directory where the t1 is
                    cfg.t1dir = [cfg.datadir,cfg.subject,'/anat/'];
                case 'Sess1_Scan2'
                    % where the unprocessed EPI 4d file is
                    cfg.rawdir = [cfg.datadir,cfg.subject,'/func_2/'];
                    % Directory where the t1 is
                    cfg.t1dir = [cfg.datadir,cfg.subject,'/anat/'];
            end

            % where the processed epi 4d files will be output to from prepro_base
            cfg.preprodir = [cfg.rawdir,'prepro/'];

            % file name of EPI 4d file
            cfg.EPI = [cfg.subject,'_task-rest_bold.nii'];
            % name of t1 file.
            cfg.t1name = [cfg.subject,'_T1w.nii'];

            % the path and filename of the template in MNI space to which everything
            % will be normalized
            cfg.mni_template = [cfg.spmdir,'templates/T1.nii'];

            % preprocessing settings
            % length of time series (no. vols)
            cfg.N = 193;
            % Repetition time of acquistion in secs
            cfg.TR = 2.46;

            % Number of slices in EPI volumes.
            cfg.numSlices = 44;
            % Vector defining acquisition order of EPI slices (necessary for slice-timing correction.
            % See help of slicetime_epis.m for guidance on how to define)
            cfg.order = [1:2:cfg.numSlices,2:2:cfg.numSlices]; % interleaved
            % Reference slice for slice timing acquisition. See SlicetimeEPI.m for guidance on how to define.
            cfg.refSlice = cfg.numSlices-1; % use for interleaved order

            % Scalar value indicating spatial smoothing kernal size in mm.
            cfg.kernel = 6;
            % Low-pass cut-off for bandpass filter in Hz (e.g., .08)
            cfg.LowPass = 10^6;
            % Hi-pass cut-off for bandpass filter in Hz (e.g., .008)
            cfg.HighPass = 0.008;
        case 'EAT'
            %% EAT project, edited from the GoC project
            % Where the subjects' directories are
            cfg.datadir = '/projects/kg98/Meadhbh/EAT_project/datadir/rawdata/';

            switch cfg.WhichSessScan
                case 'Sess1_Scan1'
                    % where the unprocessed EPI 4d file is
                    cfg.rawdir = [cfg.datadir,cfg.subject,'/func/']; %%have renamed spaces here with underscores
                    % Directory where the t1 is
                    cfg.t1dir = [cfg.datadir,cfg.subject,'/anat/'];

                    % where the processed epi 4d files will be output to from prepro_base
                    cfg.preprodir = [cfg.rawdir,'prepro_run-01/'];
                    % file name of EPI 4d file
                    cfg.EPI = [cfg.subject,'_task-eat_run-01_bold.nii'];
                case 'Sess1_Scan2'
                    % where the unprocessed EPI 4d file is
                    cfg.rawdir = [cfg.datadir,cfg.subject,'/func/']; %%have renamed spaces here with underscores
                    % Directory where the t1 is
                    cfg.t1dir = [cfg.datadir,cfg.subject,'/anat/'];

                    % where the processed epi 4d files will be output to from prepro_base
                    cfg.preprodir = [cfg.rawdir,'prepro_run-02/'];
                    % file name of EPI 4d file
                    cfg.EPI = [cfg.subject,'_task-eat_run-02_bold.nii'];
                case 'Sess1_Scan3'
                    % where the unprocessed EPI 4d file is
                    cfg.rawdir = [cfg.datadir,cfg.subject,'/func/']; %%have renamed spaces here with underscores
                    % Directory where the t1 is
                    cfg.t1dir = [cfg.datadir,cfg.subject,'/anat/'];

                    % where the processed epi 4d files will be output to from prepro_base
                    cfg.preprodir = [cfg.rawdir,'prepro_run-03/'];
                    % file name of EPI 4d file
                    cfg.EPI = [cfg.subject,'_task-eat_run-03_bold.nii'];
                case 'Sess1_Scan4'
                    % where the unprocessed EPI 4d file is
                    cfg.rawdir = [cfg.datadir,cfg.subject,'/func/']; %%have renamed spaces here with underscores
                    % Directory where the t1 is
                    cfg.t1dir = [cfg.datadir,cfg.subject,'/anat/'];

                    % where the processed epi 4d files will be output to from prepro_base
                    cfg.preprodir = [cfg.rawdir,'prepro_run-04/'];
                    % file name of EPI 4d file
                    cfg.EPI = [cfg.subject,'_task-eat_run-04_bold.nii'];
                case 'Sess1_Scan5'
                    % where the unprocessed EPI 4d file is
                    cfg.rawdir = [cfg.datadir,cfg.subject,'/func/']; %%have renamed spaces here with underscores
                    % Directory where the t1 is
                    cfg.t1dir = [cfg.datadir,cfg.subject,'/anat/'];

                    % where the processed epi 4d files will be output to from prepro_base
                    cfg.preprodir = [cfg.rawdir,'prepro_run-05/'];
                    % file name of EPI 4d file
                    cfg.EPI = [cfg.subject,'_task-eat_run-05_bold.nii'];
                case 'Sess1_Scan6'
                    % where the unprocessed EPI 4d file is
                    cfg.rawdir = [cfg.datadir,cfg.subject,'/func/']; %%have renamed spaces here with underscores
                    % Directory where the t1 is
                    cfg.t1dir = [cfg.datadir,cfg.subject,'/anat/'];

                    % where the processed epi 4d files will be output to from prepro_base
                    cfg.preprodir = [cfg.rawdir,'prepro_run-06/'];
                    % file name of EPI 4d file
                    cfg.EPI = [cfg.subject,'_task-eat_run-06_bold.nii'];
            end

            % Directory where the t1 is
            cfg.t1dir = [cfg.datadir,cfg.subject,'/anat/'];
            % name of t1 file.
            cfg.t1name = [cfg.subject,'_T1w.nii'];

            % the path and filename of the template in MNI space to which everything
            % will be normalized
            cfg.mni_template = [cfg.spmdir,'templates/T1.nii'];

            % preprocessing settings
            % length of time series (no. vols)
            cfg.N = 116; %116mb
            % Repetition time of acquistion in secs
            cfg.TR = 2.46; %must check this EAT %2460ms-mb

            % Number of slices in EPI volumes.
            cfg.numSlices = 44; %not sure why this was empty %44mb
            % Vector defining acquisition order of EPI slices (necessary for slice-timing correction.
            % See help of slicetime_epis.m for guidance on how to define)
            cfg.order = [1:2:cfg.numSlices,2:2:cfg.numSlices]; % interleaved
            % Reference slice for slice timing acquisition. See SlicetimeEPI.m for guidance on how to define.
            cfg.refSlice = cfg.numSlices-1; % use for interleaved order

            % Scalar value indicating spatial smoothing kernal size in mm.
            cfg.kernel = 6;
            % Low-pass cut-off for bandpass filter in Hz (e.g., .08)
            cfg.LowPass = 0.08;
            % Hi-pass cut-off for bandpass filter in Hz (e.g., .008)
            cfg.HighPass = 0.008;
        case 'OCDPG_DCM'
            % Where the subjects' directories are
            cfg.datadir = [cfg.parentdir_scratch,'ResProjects/rfMRI_DCM/OCDPG/data/'];

            switch cfg.WhichSessScan
                case 'Sess1_Scan1'
                    % where the unprocessed EPI 4d file is
                    cfg.rawdir = [cfg.datadir,cfg.subject,'/func/'];
                    % Directory where the t1 is
                    cfg.t1dir = [cfg.datadir,cfg.subject,'/anat/'];
            end

            % where the processed epi 4d files will be output to from prepro_base
            cfg.preprodir = [cfg.rawdir,'prepro/'];

            % file name of EPI 4d file
            cfg.EPI = [cfg.subject,'_task-rest_bold.nii'];
            % name of t1 file.
            cfg.t1name = [cfg.subject,'_T1w.nii'];

            % the path and filename of the template in MNI space to which everything
            % will be normalized
            cfg.mni_template = [cfg.spmdir,'templates/T1.nii'];

            % preprocessing settings
            % length of time series (no. vols)
            cfg.N = 189;
            % Repetition time of acquistion in secs
            cfg.TR = 2.5;
            % Number of slices in EPI volumes.
            cfg.numSlices = 44;
            % Vector defining acquisition order of EPI slices (necessary for slice-timing correction.
            % See help of slicetime_epis.m for guidance on how to define)
            cfg.order = [1:2:cfg.numSlices,2:2:cfg.numSlices]; % interleaved
            % Reference slice for slice timing acquisition. See SlicetimeEPI.m for guidance on how to define.
            cfg.refSlice = cfg.numSlices-1; % use for interleaved order

            % Scalar value indicating spatial smoothing kernal size in mm.
            cfg.kernel = 6;
            % Low-pass cut-off for bandpass filter in Hz (e.g., .08)
            cfg.LowPass = 0.08;
            % Hi-pass cut-off for bandpass filter in Hz (e.g., .008)
            cfg.HighPass = 0.008;
    end

    % ------------------------------------------------------------------------------
    % conc_TMS_FMRI
    % ------------------------------------------------------------------------------
    if ismember('conc_TMS_FMRI',cfg.WhichProject,'rows')
        if exist(cfg.t1dir) == 0
            fprintf(1,'\t\t Initialising t1dir\n')
            mkdir(cfg.t1dir)
        elseif exist(cfg.t1dir) == 7
            fprintf(1,'\t\t Cleaning and re-initialising t1dir\n')
            rmdir(cfg.t1dir,'s')
            mkdir(cfg.t1dir)
        end

        fprintf(1,'\t\t Copying t1\n')
        copyfile([cfg.datadir,cfg.subject,'/anat/',cfg.t1name,'*'],cfg.t1dir)

        if ismember('FDown_Scan1',cfg.WhichSessScan,'rows') | ismember('FDown_Scan2',cfg.WhichSessScan,'rows')
            fprintf(1,'\t\t Copying face down t1 as well...\n')
            copyfile([cfg.datadir,cfg.subject,'/anat/',cfg.t14norm,'*'],cfg.t1dir)
        end
    end

    % ------------------------------------------------------------------------------
    % run prepro_base
    runBase = 1;
    % ------------------------------------------------------------------------------
    if runBase == 1
        [cfg.tN,cfg.gm,cfg.wm,cfg.csf,cfg.epiBrainMask,cfg.t1BrainMask,cfg.BrainMask,cfg.gmmask,cfg.wmmask,cfg.csfmask,cfg.dvars,cfg.dvarsExtract,cfg.fdThr,cfg.dvarsThr,cfg.exclude,cfg.outEPI] = prepro_base(cfg);
    elseif runBase == 0
        % assumes 6P has been run
        temp = load([cfg.preprodir,'6P/cfg.mat']);
        cfg.tN = temp.cfg.tN;
        cfg.gm = temp.cfg.gm;
        cfg.wm = temp.cfg.wm;
        cfg.csf = temp.cfg.csf;
        cfg.epiBrainMask = temp.cfg.epiBrainMask;
        cfg.t1BrainMask = temp.cfg.t1BrainMask;
        cfg.BrainMask = temp.cfg.BrainMask;
        cfg.gmmask = temp.cfg.gmmask;
        cfg.wmmask = temp.cfg.wmmask;
        cfg.csfmask = temp.cfg.csfmask;
        cfg.dvars = temp.cfg.dvars;
        cfg.dvarsExtract = temp.cfg.dvarsExtract;
        cfg.fdThr = temp.cfg.fdThr;
        cfg.dvarsThr = temp.cfg.dvarsThr;
        cfg.exclude = temp.cfg.exclude;
        cfg.outEPI = temp.cfg.outEPI;
    end

    % ------------------------------------------------------------------------------
    % noise correction and time series
    % ------------------------------------------------------------------------------
    % Enter the names of the noise corrections strategies you want to use into the cell 'noiseOptions'
    % This script will then loop over each strategy
    % If you only want to run one then just use something like: noiseOptions = {'24P+aCC'};

    % All pipelines
    % noiseOptions = {'6P','6P+2P','6P+2P+GSR','24P','24P+8P','24P+8P+4GSR','12P+aCC','24P+aCC','12P+aCC50','24P+aCC50','24P+aCC+4GSR','24P+aCC50+4GSR','ICA-AROMA+2P','ICA-AROMA+2P+GSR','ICA-AROMA+8P','ICA-AROMA+8P+4GSR','24P+8P+4GSR+SpikeReg'};
    noiseOptions = {'ICA-AROMA+2P','ICA-AROMA+2P+GSR'};

    % If subject was not marked for exclusion for scrubbing, then append the JP14 pipelines
    if cfg.exclude == 0 & cfg.intnorm == 1 & cfg.runBandpass == 1
        fprintf(1, '\n\t\t Adding JP14 pipelines \n\n');
        noiseOptions2Append = {'24P+4P+2GSR+JP14Scrub'};
        noiseOptions = [noiseOptions,noiseOptions2Append];
    end

    % Loop over noise correction options
    for i = 1:length(noiseOptions)

        % Set noise correction
        cfg.removeNoise = noiseOptions{i};

        % Override smoothing order if using ICA-AROMA
        if any(~cellfun('isempty',strfind({cfg.removeNoise},'ICA-AROMA'))) == 1
            fprintf(1, '\n\t\t !!!! Forcing smoothing order to ''before'' for ICA-AROMA !!!! \n\n');
            cfg.smoothing = 'before';
        else
            cfg.smoothing = smoothing;
        end

        % define inputs to noise correction
        switch cfg.smoothing
            case {'after','none'}
                if any(~cellfun('isempty',strfind({cfg.removeNoise},'JP14Scrub'))) == 0
                    % If we run smoothing AFTER noise correction, then the input file is unsmoothed epi from prepro_base
                    cfg.CleanIn = cfg.outEPI{1};
                    % and, the nuisance inputs are the same image
                    cfg.NuisanceIn_wm = cfg.outEPI{1};
                    cfg.NuisanceIn_csf = cfg.outEPI{1};
                elseif any(~cellfun('isempty',strfind({cfg.removeNoise},'JP14Scrub'))) == 1
                    % If Power's scrubbing then its the power equivalent instead
                    cfg.CleanIn = cfg.outEPI{3};
                    cfg.NuisanceIn_wm = cfg.outEPI{3};
                    cfg.NuisanceIn_csf = cfg.outEPI{3};
                end
            case 'before'
                if any(~cellfun('isempty',strfind({cfg.removeNoise},'JP14Scrub'))) == 0
                    % If we run smoothing BEFORE noise correction, then the input file is smoothed epi from prepro_base
                    cfg.CleanIn = cfg.outEPI{2};
                    % and, the nuisance inputs are the same image
                    cfg.NuisanceIn_wm = cfg.outEPI{2};
                    cfg.NuisanceIn_csf = cfg.outEPI{2};
                elseif any(~cellfun('isempty',strfind({cfg.removeNoise},'JP14Scrub'))) == 1
                    % If Power's scrubbing then its the power equivalent instead
                    cfg.CleanIn = cfg.outEPI{4};
                    cfg.NuisanceIn_wm = cfg.outEPI{4};
                    cfg.NuisanceIn_csf = cfg.outEPI{4};
                end
        end

        % ------------------------------------------------------------------------------
        % run prepro_noise
        % ------------------------------------------------------------------------------
        [cfg.noiseTS,cfg.outdir,cfg.noiseTSz] = prepro_noise(cfg);

        % ------------------------------------------------------------------------------
        % extract time series
        runTS = 1;
        % ------------------------------------------------------------------------------
        if runTS == 1
            cd(cfg.outdir)

            cfg.parcFiles = {[cfg.parentdir,'ROIs/Gordon/Parcels_MNI_222.nii'],...
                            [cfg.parentdir,'ROIs/Power/Power.nii'],...
                            [cfg.parentdir,'ROIs/TriStri/TriStri.nii'],...
                            [cfg.parentdir,'ROIs/DiMartino/SphereParc02.nii'],...
                            [cfg.parentdir,'ROIs/HCP/MMP_in_MNI_asym_222_continuousLabels.nii'],...
                            [cfg.parentdir,'ROIs/AAL/AAL_2mm/raal.nii']};

            cfg.parcWeightGM = {'yes',...
                                'yes',...
                                'no',...
                                'no',...
                                'yes',...
                                'yes'};

            % Set input image for time series extraction
            cfg.ExtractIn = 'epi_prepro.nii';

            % Initialise roi time series variable
            cfg.roiTS = [];

            % Loop over parcellation files
            for i = 1:length(cfg.parcFiles)
                % Set parcellation file
                cfg.parcFile = cfg.parcFiles{i};
                % set GM weight
                cfg.weightGM = cfg.parcWeightGM{i};
                % extract time series
                cfg.roiTS{i} = prepro_extractTS_FSL(cfg);
            end

            cfg = rmfield(cfg,{'parcFile','weightGM'});
        end

        % Save data
        save('cfg.mat','cfg')

        % ------------------------------------------------------------------------------
        % Compress final outputs
        % ------------------------------------------------------------------------------
        fprintf('\n\t\t ----- Compressing %s outputs ----- \n\n',cfg.removeNoise);
        if any(size(dir([cfg.outdir '*.nii' ]),1))
            cd(cfg.outdir)
            gzip('*.nii')
            pause(5)
            delete('*.nii')
        end
    end

    % ------------------------------------------------------------------------------
    % Compress base & t1 outputs
    % ------------------------------------------------------------------------------
    fprintf('\n\t\t ----- Compressing base outputs ----- \n\n');
    if any(size(dir([cfg.preprodir '*.nii' ]),1)) == 1
        cd(cfg.preprodir)
        gzip('*.nii')
        pause(5)
        delete('*.nii')
    end

    if any(size(dir([cfg.t1dir '*.nii' ]),1)) == 1
        cd(cfg.t1dir)
        gzip('*.nii')
        pause(5)
        delete('*.nii')
    end

    if any(size(dir([cfg.rawdir '*.nii' ]),1)) == 1
        cd(cfg.rawdir)
        gzip('*.nii')
        pause(5)
        delete('*.nii')
    end
    fprintf('\n\t\t ----- Finished. ----- \n\n');
end
