%% HCP_run_prepro: 
function [] = HCP_run_prepro(subject)
    cfg = struct;
    cfg.subject = subject;

    % ------------------------------------------------------------------------------
    % Store date and time
    % ------------------------------------------------------------------------------
    cfg.DateTime = datetime('now');
    
    % ------------------------------------------------------------------------------
    % Add paths - edit this section
    % ------------------------------------------------------------------------------
        % where the prepro scripts are
        cfg.scriptdir = '/home/lindenmp/kg98/Linden/Scripts/rs-fMRI/prepro/';
        addpath(cfg.scriptdir)
        cfg.funcdir = '/home/lindenmp/kg98/Linden/Scripts/rs-fMRI/func/';
        addpath(cfg.funcdir)

        % where spm is
        cfg.spmdir = '/usr/local/spm8/matlab2015b.r6685/';
        addpath(cfg.spmdir)

        % set FSL environments 
        cfg.fsldir = '/usr/local/fsl/5.0.9/fsl/bin/';
        setenv('FSLDIR',cfg.fsldir(1:end-4));
        setenv('FSLOUTPUTTYPE','NIFTI');
        setenv('LD_LIBRARY_PATH',[getenv('PATH'),getenv('LD_LIBRARY_PATH'),':/usr/lib/fsl/5.0'])

    % ------------------------------------------------------------------------------
    % Set project settings and parameters
    % Use WhichProject if you're juggling multiple datasets
    % ------------------------------------------------------------------------------
        cfg.projdir = '/projects/kg98/Linden/ResProjects/HCP_BF_TimeSeries/';
        % Where the subjects' directories are
        cfg.datadir = [cfg.projdir,'data/'];

        % Setup output directory
        cfg.preprodir = ([cfg.datadir,cfg.subject,'/']); 
        if exist(cfg.preprodir) == 0
            fprintf(1,'Initialising cfg.preprodir\n')
            mkdir(cfg.preprodir)
        elseif exist(cfg.preprodir) == 7
            fprintf(1,'Cleaning and re-initialising cfg.preprodir\n')
            rmdir(cfg.preprodir,'s')
            mkdir(cfg.preprodir)
        end

        cfg.HCPdir = '/scratch/hcp/';

        % where the unprocessed EPI 4d file is
        cfg.rawdir = [cfg.HCPdir,cfg.subject,'/MNINonLinear/Results/rfMRI_REST1_LR/'];

        % file name of EPI 4d file
        cfg.EPI = 'rfMRI_REST1_LR_hp2000_clean.nii.gz';

        cfg.roidir = [cfg.HCPdir,cfg.subject,'/MNINonLinear/ROIs/'];
        cfg.parcFile = 'wmparc.2.nii.gz';

        % Directory where the t1 is
        cfg.t1dir = [cfg.HCPdir,cfg.subject,'/MNINonLinear/']; 
        % name of t1 file.
        cfg.t1name = 'T1w.nii.gz';

        % preprocessing settings
        % length of time series (no. vols)
        cfg.N = 1200; 
        % Repetition time of acquistion in secs
        cfg.TR = 0.720;
        % Scalar value indicating spatial smoothing kernal size in mm. 
        % cfg.kernel = 8;
        % Low-pass cut-off for bandpass filter in Hz (e.g., .08) 
        cfg.LowPass = 10^6;
        % Hi-pass cut-off for bandpass filter in Hz (e.g., .008)
        cfg.HighPass = 0.008;

    % ------------------------------------------------------------------------------
    % Copy input files over to kg98
    % ------------------------------------------------------------------------------
        % T1
        copyfile([cfg.t1dir,cfg.t1name],cfg.preprodir)
        % EPI
        copyfile([cfg.rawdir,cfg.EPI],cfg.preprodir)
        % mov params
        copyfile([cfg.rawdir,'/Movement*'],cfg.preprodir)


    % ------------------------------------------------------------------------------
    % Unzip
    % ------------------------------------------------------------------------------
        cd(cfg.preprodir)
        % T1
        gunzip(cfg.t1name); delete(cfg.t1name); cfg.t1name = cfg.t1name(1:end-3);
        % EPI
        gunzip(cfg.EPI); delete(cfg.EPI); cfg.EPI = cfg.EPI(1:end-3);

    % ------------------------------------------------------------------------------
    % Segment T1
    % ------------------------------------------------------------------------------
        cd(cfg.preprodir)
        % Tissue segment T1 with SPM
        SegmentT1([cfg.preprodir,cfg.t1name],cfg.spmdir,0,0);

        % outputs
        gm = ['c1',cfg.t1name];
        wm = ['c2',cfg.t1name];
        csf = ['c3',cfg.t1name];

        % reslice to EPI dimensions
        imat = eye(4); dlmwrite('eye.mat',imat,'\t')
        system([cfg.fsldir,'flirt -in ',gm,' -ref ',cfg.EPI,' -applyxfm -init eye.mat -out ',gm]);
        system([cfg.fsldir,'flirt -in ',wm,' -ref ',cfg.EPI,' -applyxfm -init eye.mat -out ',wm]);
        system([cfg.fsldir,'flirt -in ',csf,' -ref ',cfg.EPI,' -applyxfm -init eye.mat -out ',csf]);

    % ------------------------------------------------------------------------------
    % Create binary brain mask
    % ------------------------------------------------------------------------------
        fprintf(1,'\n\t\t ----- Create brain masks ----- \n\n')

        % EPI
        cd(cfg.preprodir)
        MaskIn = cfg.EPI;
        system([cfg.fsldir,'bet ',MaskIn,' epi_brain -f 0.4 -n -m -R']);
        delete('epi_brain.nii')
        cfg.epiBrainMask = 'epi_brain_mask.nii';

        % T1
        MaskIn = cfg.t1name;
        system([cfg.fsldir,'bet ',MaskIn,' t1_brain -f 0.4 -n -m -R']);
        delete('t1_brain.nii')
        cfg.t1BrainMask = 't1_brain_mask.nii';
        % reslice
        system([cfg.fsldir,'flirt -in ',cfg.t1BrainMask,' -ref ',cfg.EPI,' -interp nearestneighbour -applyxfm -init eye.mat -out ',cfg.t1BrainMask]);

        % Union
        system([cfg.fsldir,'fslmaths ',cfg.epiBrainMask,' -add ',cfg.t1BrainMask,' -bin brain_mask']);
        cfg.BrainMask = 'brain_mask.nii';

    % ------------------------------------------------------------------------------
    % Mask out non-brain tissue from EPI image, T1, and tissue masks
    % ------------------------------------------------------------------------------
        fprintf(1,'\n\t\t ----- Apply brain masks ----- \n\n')
        
        cd(cfg.preprodir)
        
        % Tissue masks
        system([cfg.fsldir,'fslmaths ',gm,' -mas ',cfg.preprodir,cfg.BrainMask,' b',gm]);
        system([cfg.fsldir,'fslmaths ',wm,' -mas ',cfg.preprodir,cfg.BrainMask,' b',wm]);
        system([cfg.fsldir,'fslmaths ',csf,' -mas ',cfg.preprodir,cfg.BrainMask,' b',csf]);

        % T1 outputs
        gm = ['b',gm];
        wm = ['b',wm];
        csf = ['b',csf];

    % ------------------------------------------------------------------------------
    % Generate MNI tissue masks for noise correction and tissue smoothing
    % ------------------------------------------------------------------------------
        cd(cfg.preprodir)

        % retain only voxels with >=50% probability of gm
        % Use this for gm blurring
        system([cfg.fsldir,'fslmaths ',gm,' -thr .50 -bin gm50_bin']);

        % retain anything with >=1% probability of being gm
        % This is just to remove gm overlap in wm/csf masks
        system([cfg.fsldir,'fslmaths ',gm,' -thr .01 -bin gm01_bin']);
        % retain only voxels with >=99% probability of wm
        system([cfg.fsldir,'fslmaths ',wm,' -thr .99 -bin wm99_bin']);
        % retain only voxels with >=99% probability of csf
        system([cfg.fsldir,'fslmaths ',csf,' -thr .99 -bin csf99_bin']);

        % remove overlap between gm and white and csf masks
        system([cfg.fsldir,'fslmaths gm01_bin -mul -1 -add 1 gm01_inv']);
        system([cfg.fsldir,'fslmaths wm99_bin -mul gm01_inv wm_final']);
        system([cfg.fsldir,'fslmaths csf99_bin -mul gm01_inv csf_final']);

        % erode
        % system([cfg.fsldir,'fslmaths wm_final -ero wm_final']);
        % system([cfg.fsldir,'fslmaths csf_final -ero csf_final']);

        cfg.gmmask = 'gm50_bin.nii';
        cfg.wmmask = 'wm_final.nii';
        cfg.csfmask = 'csf_final.nii';

    % ------------------------------------------------------------------------------
    % Detrend EPI
    % ------------------------------------------------------------------------------
        fprintf(1,'\n\t\t ----- Detrending ----- \n\n')

        DetrendIn = cfg.EPI;
        DetrendOut = ['d',DetrendIn];

        cd(cfg.preprodir)
        % create separate directory for REST detrend function
        dtdir = [cfg.preprodir,'temp'];
        mkdir(dtdir)
        % move 4d file to directory
        movefile([cfg.preprodir,DetrendIn],dtdir,'f')
        
        % Split 4D file in 3D
        cd(dtdir)
        % spm_file_split(DetrendIn)
        system([cfg.fsldir,'fslsplit ',cfg.EPI]);

        % Move 4D file back to cfg.preprodir
        movefile(DetrendIn,cfg.preprodir,'f')

        % detrend epis using rest
        cd(cfg.preprodir)
        rest_detrend(dtdir, '_detrend')

        % delete 3D files
        cd(dtdir)
        delete('*nii')

        % move output file back to cfg.preprodir
        movefile([dtdir,'_detrend/detrend_4DVolume.nii'],[cfg.preprodir,DetrendOut],'f')

        cd(cfg.preprodir)
        rmdir('temp','s')
        rmdir('temp_detrend','s')

    % ------------------------------------------------------------------------------
    % Noise correction
    % ------------------------------------------------------------------------------
        cfg.CleanIn = DetrendOut;
        cfg.removeNoise = '24P+2P+GSR';

        % Setup preprodir outdir
        cfg.outdir = [cfg.preprodir,cfg.removeNoise,'/']; 
        if exist(cfg.outdir) == 0
            mkdir(cfg.outdir)
        end

    % ------------------------------------------------------------------------------
    % 1) Motion parameters
    % ------------------------------------------------------------------------------
        % read in motion
        mfile = dir([cfg.preprodir,'Movement_Regressors_dt.txt']);
        mov = dlmread([cfg.preprodir,mfile(1).name]);

    % ------------------------------------------------------------------------------
    % 2) Physiological time series
    % ------------------------------------------------------------------------------
        cd(cfg.outdir)
        % Extract time series

        % White Matter
        str = [cfg.fsldir,'fslmeants -i ',cfg.preprodir,cfg.CleanIn,' -o wmTS.txt -m ',cfg.preprodir,cfg.wmmask];
        system(str); 
        wmTS = dlmread('wmTS.txt');

        % CSF
        str = [cfg.fsldir,'fslmeants -i ',cfg.preprodir,cfg.CleanIn,' -o csfTS.txt -m ',cfg.preprodir,cfg.csfmask];
        system(str); 
        csfTS = dlmread('csfTS.txt');

    % ------------------------------------------------------------------------------
    % 2.1) Global Signal Regression
    % ------------------------------------------------------------------------------
        % Global
        str = [cfg.fsldir,'fslmeants -i ',cfg.preprodir,cfg.CleanIn,' -o gsTS.txt -m ',cfg.preprodir,cfg.BrainMask];
        system(str); 
        gsTS = dlmread('gsTS.txt');

    % ------------------------------------------------------------------------------
    % Generate noiseTS variable for fsl_regfilt
    % ------------------------------------------------------------------------------
        fprintf(1,'\n\t\t ----- Generating noiseTS ----- \n\n');

        noiseTS = mov(:,1:6);
        noiseTS = GetDerivatives(noiseTS);
        noiseTS = [noiseTS wmTS csfTS gsTS];

    % ------------------------------------------------------------------------------
    % Clean data: fsl_regfilt
    % ------------------------------------------------------------------------------
        cd(cfg.outdir)
        fprintf(1,'\n\t\t ----- Running nuisance regression ----- \n\n');
        
        % also, write out noise signals as text file (for regfilt)
        dlmwrite('noiseTS.txt',noiseTS,'delimiter','\t','precision','%.6f');

        % clean data with fsl_regfilt
        % Linden: x is a variable that stores the -f flag input for fsl_regfilt.
        %     e.g., -f 1,2,3,4,5
        x = regexprep(num2str(1:size(noiseTS,2)),' ',',');
        x = regexprep(x,',,,',',');
        x = regexprep(x,',,',',');
        x = ['"',x,'"'];

        str = [cfg.fsldir,'fsl_regfilt -i ',cfg.preprodir,cfg.CleanIn,' -o epi_clean -d noiseTS.txt -f ',x];
        system(str);

        clear x

        system('gunzip -rf *epi_clean*');

    % ------------------------------------------------------------------------------
    % Bandpass filter with REST
    % ------------------------------------------------------------------------------
        fprintf(1,'\n\t\t ----- Running bandpass filtering ----- \n\n');

        cfg.FiltIn = 'epi_clean.nii';
        CUTNUMBER = 1;

        cd(cfg.outdir)
        % creater input directory for REST function
        mkdir('temp');
        filtdir = [cfg.outdir,'temp/'];
        movefile(cfg.FiltIn,filtdir,'f');

        % split EPI into 3D files
        cd(filtdir)
        % spm_file_split(cfg.FiltIn)
        system([cfg.fsldir,'fslsplit ',cfg.FiltIn]);

        % Move 4D file back to cfg.outdir
        movefile(cfg.FiltIn,cfg.outdir,'f')

        cd(cfg.outdir)
        % bandpass epis using rest
        rest_bandpass(filtdir,cfg.TR,cfg.LowPass,cfg.HighPass,'No',[cfg.preprodir,cfg.BrainMask],CUTNUMBER)

        % Move output files back to cfg.outdir
        rmdir(filtdir,'s')
        movefile('temp_filtered/*',cfg.outdir)
        rmdir('temp_filtered','s')

        % Rename
        movefile('Filtered_4DVolume.nii','epi_prepro.nii')

    % ------------------------------------------------------------------------------
    % get time series
    % ------------------------------------------------------------------------------
        idx = dlmread([cfg.projdir,'idx2retain.txt']);
        cfg.parcFile = [cfg.HCPdir,cfg.subject,'/MNINonLinear/ROIs/wmparc.2.nii.gz'];
        cfg.ExtractIn = [cfg.preprodir,cfg.removeNoise,'/epi_prepro.nii'];
        cfg.weightGM = 'no';

        cfg.roiTS{1} = prepro_extractTS_FSL(cfg);
        cfg.roiTS{1} = cfg.roiTS{1}(:,idx);

        cd(cfg.outdir)
        save('cfg.mat','cfg')
end
