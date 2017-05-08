%% run_prepro: 
function [] = run_prepro(WhichProject,WhichSessScan,subject,discard,slicetime,despike,smoothing)
    cfg.WhichProject = WhichProject;
    cfg.WhichSessScan = WhichSessScan;
    cfg.subject = subject;
    % preprocessing options
    cfg.discard = discard;
    cfg.slicetime = slicetime;
    cfg.despike = despike;
    cfg.smoothing = smoothing;

    % ------------------------------------------------------------------------------
    % Add paths - edit this section
    % ------------------------------------------------------------------------------
        % where the prepro scripts are
        cfg.scriptdir = '/gpfs/M2Home/projects/Monash076/Linden/scripts/rfMRI-PrePro/';
        addpath(cfg.scriptdir)
        cfg.funcdir = '/gpfs/M2Home/projects/Monash076/Linden/scripts/rfMRI-Func/';
        addpath(cfg.funcdir)

        % where spm is
        cfg.spmdir = '/usr/local/spm8/matlab2014a.r5236/';
        addpath(cfg.spmdir)

        % set FSL environments 
        cfg.fsldir = '/usr/local/fsl/5.0.9/bin/';
        setenv('FSLDIR',cfg.fsldir(1:end-4));
        setenv('FSLOUTPUTTYPE','NIFTI');
        setenv('LD_LIBRARY_PATH',[getenv('PATH'),getenv('LD_LIBRARY_PATH'),':/usr/lib/fsl/5.0'])

        % ANTs
        cfg.antsdir = '/usr/local/ants/1.9.v4/bin/';
        setenv('ANTSPATH',cfg.antsdir);

        % Directory to AFNI functions
        cfg.afnidir = '/usr/local/afni/16.2.16/';

    % ------------------------------------------------------------------------------
    % Set project settings and parameters
    % Use WhichProject if you're juggling multiple datasets
    % ------------------------------------------------------------------------------
    switch cfg.WhichProject
        case 'OCDPG'
            % Where the subjects' directories are
            cfg.datadir = '/gpfs/M2Home/projects/Monash076/Linden/OCDPG/data/';

            switch cfg.WhichSessScan
                case 'Sess1_Scan1'
                    % where the unprocessed EPI 4d file is
                    cfg.rawdir = [cfg.datadir,cfg.subject,'/rfMRI/'];
                    % Directory where the t1 is
                    cfg.t1dir = [cfg.datadir,cfg.subject,'/t1/']; 
            end
            
            % where the processed epi 4d files will be output to from prepro_base
            cfg.preprodir = [cfg.rawdir,'prepro/'];

            % file name of EPI 4d file
            cfg.EPI = 'epi.nii';
            % name of t1 file.
            cfg.t1name = 't1.nii';

            % the path and filename of the template in MNI space to which everything
            % will be normalized
            cfg.mni_template = [cfg.spmdir,'templates/T1.nii'];    

            % preprocessing settings
            % length of time series (no. vols)
            cfg.N = 189; 
            % Repetition time of acquistion in secs
            cfg.TR = 2.5;
            % Desired voxel dimension (in mm) of analysis after spatial normalization
            cfg.voxdim = 2;
            % Number of slices in EPI volumes.
            cfg.numSlices = 44;
            % Vector defining acquisition order of EPI slices (necessary for slice-timing correction.
            % See help of slicetime_epis.m for guidance on how to define)
            % cfg.order = [1:1:cfg.numSlices]; % ascending
            % cfg.order = [cfg.numSlices:-1:1]; % descending
            cfg.order = [1:2:cfg.numSlices,2:2:cfg.numSlices]; % interleaved
            % cfg.order = [2:2:cfg.numSlices, 1:2:cfg.numSlices-1]; %interleave alt
            % Reference slice for slice timing acquisition. See SlicetimeEPI.m for guidance on how to define. 
            % cfg.refSlice = round(numSlices/2)); % use for sequential order (e.g., ascending or descending)
            cfg.refSlice = cfg.numSlices-1; % use for interleaved order
            % cfg.refSlice = cfg.numSlices; % use for interleaved alt order

            % Scalar value indicating spatial smoothing kernal size in mm. 
            cfg.kernel = 8;
            % Low-pass cut-off for bandpass filter in Hz (e.g., .08) 
            cfg.LowPass = 0.08;
            % Hi-pass cut-off for bandpass filter in Hz (e.g., .008)
            cfg.HighPass = 0.008;
            
            % Select whether to run certain functions in 4D or 3D mode.
                % Note that BOTH options assume the files are stored in 4D.
                % I have found that SPM sometimes does not handle loading in 4D files 
                % that consist of lots of volumes - e.g., in multiband cases or very long sequences (400+ volumes).
                % In which case, using the '3D' option will take your 4D file, split it up, do something, then concatenate back.
                % If you've just got a fairly typical fMRI run, leave as '4D'
            cfg.WhichNii = '4D';
        case 'UCLA'
            % Where the subjects' directories are
            cfg.datadir = '/gpfs/M2Home/projects/Monash076/Linden/UCLA/data/';

            switch cfg.WhichSessScan
                case 'Sess1_Scan1'
                    % where the unprocessed EPI 4d file is
                    cfg.rawdir = [cfg.datadir,cfg.subject,'/rfMRI/'];
                    % Directory where the t1 is
                    cfg.t1dir = [cfg.datadir,cfg.subject,'/t1/']; 
            end

            % where the processed epi 4d files will be output to from prepro_base
            cfg.preprodir = [cfg.rawdir,'prepro/']; 
            
            % file name of EPI 4d file
            cfg.EPI = 'epi.nii';
            % name of t1 file.
            cfg.t1name = 't1.nii';

            % the path and filename of the template in MNI space to which everything
            % will be normalized
            cfg.mni_template = [cfg.spmdir,'templates/T1.nii'];    

            % preprocessing settings
            % length of time series (no. vols)
            cfg.N = 152; 
            % Repetition time of acquistion in secs
            cfg.TR = 2;
            % Desired voxel dimension (in mm) of analysis after spatial normalization
            cfg.voxdim = 2;
            % Number of slices in EPI volumes.
            cfg.numSlices = 34;
            % Vector defining acquisition order of EPI slices (necessary for slice-timing correction.
            % See help of slicetime_epis.m for guidance on how to define)
            % cfg.order = [1:1:cfg.numSlices]; % ascending
            % cfg.order = [cfg.numSlices:-1:1]; % descending
            cfg.order = [1:2:cfg.numSlices,2:2:cfg.numSlices]; % interleaved
            % Reference slice for slice timing acquisition. See SlicetimeEPI.m for guidance on how to define. 
            % refSlice = round(numSlices/2)); % use for sequential order (e.g., ascending or descending)
            cfg.refSlice = cfg.numSlices-1; % use for interleaved order
            % Scalar value indicating spatial smoothing kernal size in mm. 
            cfg.kernel = 8;
            % Low-pass cut-off for bandpass filter in Hz (e.g., .08) 
            cfg.LowPass = 0.08;
            % Hi-pass cut-off for bandpass filter in Hz (e.g., .008)
            cfg.HighPass = 0.008;
            
            % Select whether to run certain functions in 4D or 3D mode.
                % Note that BOTH options assume the files are stored in 4D.
                % I have found that SPM sometimes does not handle loading in 4D files 
                % that consist of lots of volumes - e.g., in multiband cases or very long sequences (400+ volumes).
                % In which case, using the '3D' option will take your 4D file, split it up, do something, then concatenate back.
                % If you've just got a fairly typical fMRI run, leave as '4D'
            cfg.WhichNii = '4D';
        case 'NYU_2'
            % Where the subjects' directories are
            cfg.datadir = '/gpfs/M2Home/projects/Monash076/Linden/NYU_2/data/';
            
            % where the unprocessed EPI 4d file is
            % cfg.WhichSessScan = 'Sess1_Scan1';
            switch cfg.WhichSessScan
                case 'Sess1_Scan1'
                    % where the unprocessed EPI 4d file is
                    cfg.rawdir = [cfg.datadir,cfg.subject,'/session_1/rest_1/'];
                    % Directory where the t1 is
                    cfg.t1dir = [cfg.datadir,cfg.subject,'/session_1/anat_1/']; 
                case 'Sess1_Scan2'
                    cfg.rawdir = [cfg.datadir,cfg.subject,'/session_1/rest_2/'];
                    cfg.t1dir = [cfg.datadir,cfg.subject,'/session_1/anat_1/']; 
                case 'Sess2_Scan1'
                    cfg.rawdir = [cfg.datadir,cfg.subject,'/session_2/rest_1/'];
                    cfg.t1dir = [cfg.datadir,cfg.subject,'/session_2/anat_1/']; 
            end

            % where the processed epi 4d files will be output to from prepro_base
            cfg.preprodir = [cfg.rawdir,'prepro/'];
            
            % file name of EPI 4d file
            cfg.EPI = 'rest.nii';
            % name of t1 file.
            cfg.t1name = 'anat.nii';

            % the path and filename of the template in MNI space to which everything
            % will be normalized
            cfg.mni_template = [cfg.spmdir,'templates/T1.nii'];    

            % preprocessing settings
            % length of time series (no. vols)
            cfg.N = 180; 
            % Repetition time of acquistion in secs
            cfg.TR = 2;
            % Desired voxel dimension (in mm) of analysis after spatial normalization
            cfg.voxdim = 2;
            % Number of slices in EPI volumes.
            cfg.numSlices = 33;
            % Vector defining acquisition order of EPI slices (necessary for slice-timing correction.
            % See help of slicetime_epis.m for guidance on how to define)
            % cfg.order = [1:1:cfg.numSlices]; % ascending
            % cfg.order = [cfg.numSlices:-1:1]; % descending
            cfg.order = [1:2:cfg.numSlices,2:2:cfg.numSlices]; % interleaved
            % cfg.order = [2:2:cfg.numSlices, 1:2:cfg.numSlices-1]; %interleave alt
            % Reference slice for slice timing acquisition. See SlicetimeEPI.m for guidance on how to define. 
            % cfg.refSlice = round(numSlices/2)); % use for sequential order (e.g., ascending or descending)
            cfg.refSlice = cfg.numSlices-1; % use for interleaved order
            % cfg.refSlice = cfg.numSlices; % use for interleaved alt order

            % Scalar value indicating spatial smoothing kernal size in mm. 
            cfg.kernel = 8;
            % Low-pass cut-off for bandpass filter in Hz (e.g., .08) 
            cfg.LowPass = 0.08;
            % Hi-pass cut-off for bandpass filter in Hz (e.g., .008)
            cfg.HighPass = 0.008;
            
            % Select whether to run certain functions in 4D or 3D mode.
                % Note that BOTH options assume the files are stored in 4D.
                % I have found that SPM sometimes does not handle loading in 4D files 
                % that consist of lots of volumes - e.g., in multiband cases or very long sequences (400+ volumes).
                % In which case, using the '3D' option will take your 4D file, split it up, do something, then concatenate back.
                % If you've just got a fairly typical fMRI run, leave as '4D'
            cfg.WhichNii = '4D';
        case 'GoC'
            % Where the subjects' directories are
            cfg.datadir = '/projects/kg98/kristina/GenofCog/data/';

             switch cfg.WhichSessScan
                case 'Sess1_Scan1'
                    % where the unprocessed EPI 4d file is
                    cfg.rawdir = [cfg.datadir,cfg.subject,'/rfMRI/'];
                    % Directory where the t1 is
                    cfg.t1dir = [cfg.datadir,cfg.subject,'/t1/']; 
            end


            % where the processed epi 4d files will be output to from prepro_base
            cfg.preprodir = [cfg.rawdir,'prepro/'];
            % file name of EPI 4d file
            cfg.EPI = 'epi.nii';

            % Directory where the t1 is
            cfg.t1dir = [cfg.datadir,cfg.subject,'/t1/']; 
            % name of t1 file.
            cfg.t1name = 't1.nii';

            % the path and filename of the template in MNI space to which everything
            % will be normalized
            cfg.mni_template = [cfg.spmdir,'templates/T1.nii'];    

            % preprocessing settings
            % length of time series (no. vols)
            cfg.N = 620; 
            % Repetition time of acquistion in secs
            cfg.TR = 0.754;
            % Desired voxel dimension (in mm) of analysis after spatial normalization
            cfg.voxdim = 2;
            % Number of slices in EPI volumes.
            cfg.numSlices = [];
            % Vector defining acquisition order of EPI slices (necessary for slice-timing correction.
            % See help of slicetime_epis.m for guidance on how to define)
            % cfg.order = [1:1:cfg.numSlices]; % ascending
            % cfg.order = [cfg.numSlices:-1:1]; % descending
            cfg.order = []; % interleaved
            % Reference slice for slice timing acquisition. See SlicetimeEPI.m for guidance on how to define. 
            % cfg.refSlice = round(numSlices/2)); % use for sequential order (e.g., ascending or descending)
            cfg.refSlice = []; % use for interleaved order
            % Scalar value indicating spatial smoothing kernal size in mm. 
            cfg.kernel = 8;
            % Low-pass cut-off for bandpass filter in Hz (e.g., .08) 
            cfg.LowPass = 0.08;
            % Hi-pass cut-off for bandpass filter in Hz (e.g., .008)
            cfg.HighPass = 0.008;

            % Select whether to run certain functions in 4D or 3D mode.
                % Note that BOTH options assume the files are stored in 4D.
                % I have found that SPM sometimes does not handle loading in 4D files 
                % that consist of lots of volumes - e.g., in multiband cases or very long sequences (400+ volumes).
                % In which case, using the '3D' option will take your 4D file, split it up, do something, then concatenate back.
                % If you've just got a fairly typical fMRI run, leave as '4D'
            cfg.WhichNii = '3D';
    end

    % ------------------------------------------------------------------------------
    % run prepro_base
    runBase = 1;
    % ------------------------------------------------------------------------------
    if runBase == 1
        [cfg.tN,cfg.gm,cfg.wm,cfg.csf,cfg.epiBrainMask,cfg.t1BrainMask,cfg.BrainMask,cfg.gmmask,cfg.wmmask,cfg.csfmask,cfg.outEPI] = prepro_base(cfg);
    elseif runBase == 0
        % Alot of this is hard coded just for testing purposes.
        % Typically users will just set runBase = 1;
        cfg.tN = cfg.N - 4;
        cfg.gm = ['bwc1c',cfg.t1name];
        cfg.wm = ['bwc2c',cfg.t1name];
        cfg.csf = ['bwc3c',cfg.t1name];
        cfg.epiBrainMask = 'epi_brain_mask.nii';
        cfg.t1BrainMask = 't1_brain_mask.nii';
        cfg.BrainMask = 'brain_mask.nii';
        cfg.gmmask = 'gm50_bin.nii';
        cfg.wmmask = 'wm_final.nii';
        cfg.csfmask = 'csf_final.nii';
        cfg.outEPI{1} = ['idbwrdat',cfg.EPI];
        cfg.outEPI{2} = ['s_idbwrdat',cfg.EPI];
        cfg.outEPI{3} = ['sgm_idbwrdat',cfg.EPI];
        cfg.outEPI{4} = ['swm_idbwrdat',cfg.EPI];
        cfg.outEPI{5} = ['scsf_idbwrdat',cfg.EPI];
    end

    % ------------------------------------------------------------------------------
    % noise correction and time series
    % ------------------------------------------------------------------------------
    % Enter the names of the noise corrections strategies you want to use into the cell 'noiseOptions'
    % This script will then loop over each strategy
    % If you only want to run one then just use something like: noiseOptions = {'24P+aCC'};
    
    % All pipelines
    noiseOptions = {'6P','6P+2P','6P+2P+GSR','24P','24P+8P','24P+8P+4GSR','24P+8P+SpikeReg','24P+8P+4GSR+SpikeReg','12P+aCC','24P+aCC','12P+aCC50','24P+aCC50','24P+aCC+4GSR','24P+aCC50+4GSR','24P+aCC+SpikeReg','24P+aCC+4GSR+SpikeReg','ICA-AROMA+2P','ICA-AROMA+2P+SpikeReg','ICA-AROMA+GSR','ICA-AROMA+2P+GSR','ICA-AROMA+8P','ICA-AROMA+4GSR','ICA-AROMA+8P+4GSR'};

    % Loop over noise correction options
    for i = 1:length(noiseOptions)

        % Incase you want to test different things, use suffix variable to manually name output directory
        % set to cfg.suffix = ''; if no suffix is needed
        cfg.suffix = '';
        
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
            case 'before'
                % If we run smoothing BEFORE noise correction, then the input file is the smoothed gm epi from prepro_base
                cfg.CleanIn = cfg.outEPI{2};
                % and, the nuisance inputs are the tissue-specific smoothed wm and csf epi images
                cfg.NuisanceIn_wm = cfg.outEPI{4};
                cfg.NuisanceIn_csf = cfg.outEPI{5};
            case {'after','none'}
                % If we run smoothing AFTER noise correction, then the input file is unsmoothed epi from prepro_base
                cfg.CleanIn = cfg.outEPI{1};
                % and, the nuisance inputs are the same image
                cfg.NuisanceIn_wm = cfg.outEPI{1};
                cfg.NuisanceIn_csf = cfg.outEPI{1};
        end

        % ------------------------------------------------------------------------------
        % run prepro_noise
        runNoise = 1;
        % ------------------------------------------------------------------------------
        if runNoise == 1
            [cfg.noiseTS,cfg.outdir] = prepro_noise(cfg);
        elseif runNoise == 0
            if cfg.CleanIn(1) == 's'
                cfg.outdir = [cfg.preprodir,'s',cfg.removeNoise,cfg.suffix,'/'];
            else
                cfg.outdir = [cfg.preprodir,cfg.removeNoise,cfg.suffix,'/'];
            end
            cfg.noiseTS = dlmread([cfg.outdir,'noiseTS.txt']);
        end

        % ------------------------------------------------------------------------------
        % extract time series
        runTS = 1;
        % ------------------------------------------------------------------------------
        if runTS == 1
            cd(cfg.outdir)
            
            % Parcellation file for time series extraction
            cfg.parcFiles = {'/gpfs/M2Home/projects/Monash076/Linden/ROIs/Gordon/Parcels_MNI_222.nii',...
                            '/gpfs/M2Home/projects/Monash076/Linden/ROIs/Power/Power.nii',...
                            '/gpfs/M2Home/projects/Monash076/Linden/ROIs/TriStri/TriStri.nii',...
                            '/gpfs/M2Home/projects/Monash076/Linden/ROIs/DiMartino/SphereParc02.nii'};

            cfg.parcWeightGM = {'yes',...
                                'yes',...
                                'no',...
                                'no'};

            % Set input image for time series extraction
            switch cfg.smoothing
                case {'before','none'}
                    % If smoothing was done before noise correction:
                    cfg.ExtractIn = 'epi_prepro.nii';
                case 'after'
                    % If it was done after:
                    cfg.ExtractIn = 'sepi_prepro.nii';
            end

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
    end

end