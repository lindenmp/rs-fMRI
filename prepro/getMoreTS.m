function [] = getMoreTS(WhichProject,subject)
    % ------------------------------------------------------------------------------
    % Add paths - edit this section
    % ------------------------------------------------------------------------------
        % where the prepro scripts are
        cfg.scriptdir = '/gpfs/M2Home/projects/Monash076/Linden/scripts/rs-fMRI/prepro/';
        addpath(cfg.scriptdir)
        cfg.funcdir = '/gpfs/M2Home/projects/Monash076/Linden/scripts/rs-fMRI/func/';
        addpath(cfg.funcdir)

        % set FSL environments 
        cfg.fsldir = '/usr/local/fsl/5.0.9/bin/';
        setenv('FSLDIR',cfg.fsldir(1:end-4));
        setenv('FSLOUTPUTTYPE','NIFTI');
        setenv('LD_LIBRARY_PATH',[getenv('PATH'),getenv('LD_LIBRARY_PATH'),':/usr/lib/fsl/5.0'])

    % ------------------------------------------------------------------------------
    % Set project settings and parameters
    % Use WhichProject if you're juggling multiple datasets
    % ------------------------------------------------------------------------------
    switch cfg.WhichProject
        case 'OCDPG'
            % Where the subjects' directories are
            cfg.datadir = '/gpfs/M2Home/projects/Monash076/Linden/OCDPG/data/';

            % where the processed epi 4d files will be output to from prepro_base
            cfg.preprodir = [cfg.rawdir,'prepro/'];
    end

    
    % All pipelines
    noiseOptions = {'6P','6P+2P','6P+2P+GSR','24P','24P+8P','24P+8P+4GSR','24P+8P+SpikeReg','24P+8P+4GSR+SpikeReg','12P+aCC','24P+aCC','12P+aCC50','24P+aCC50','24P+aCC+4GSR','24P+aCC50+4GSR','24P+aCC+SpikeReg','24P+aCC+4GSR+SpikeReg','ICA-AROMA+2P','ICA-AROMA+2P+SpikeReg','ICA-AROMA+GSR','ICA-AROMA+2P+GSR','ICA-AROMA+8P','ICA-AROMA+4GSR','ICA-AROMA+8P+4GSR'};

    % Loop over noise correction options
    for i = 1:length(noiseOptions)
        removeNoise = noiseOptions{i};

        cd([cfg.datadir,subject,'/rfMRI/prepro/',removeNoise])

        load('cfg.mat')

        % ------------------------------------------------------------------------------
        % extract time series
        runTS = 1;
        % ------------------------------------------------------------------------------
        if runTS == 1
            cd(cfg.outdir)
            
            % Parcellation file for time series extraction
            cfg.new_parcFiles = {'/gpfs/M2Home/projects/Monash076/Linden/ROIs/TriStri/Striatum.nii'};

            cfg.new_parcWeightGM = {'no'};

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
            % cfg.roiTS = [];

            % Loop over parcellation files
            for i = 1:length(cfg.new_parcFiles)
                % Set parcellation file
                cfg.parcFile = cfg.new_parcFiles{i};
                % set GM weight
                cfg.weightGM = cfg.new_parcWeightGM{i};
                % extract time series 
                cfg.roiTS{end+1} = prepro_extractTS_FSL(cfg);
            end

            cfg = rmfield(cfg,{'parcFile','weightGM'});
        end

        % Save data
        save('cfg.mat','cfg')
    end

end