function [noiseTS,outdir] = prepro_noise(cfg)
    % 1 - Correct noise in EPI output from step 11 of prepro_base script (or output from step 12 in case of ICA-AROMA)
    %       See below for details
    % 2 - Bandpass filter (includes demeaning)
    % 3 - Spatial smoothing (not in case of ICA-AROMA)
    % 
    % Linden Parkes, Brain & Mental Health Laboratory, 2016
    % ------------------------------------------------------------------------------
    % Choose noise removal
    % ------------------------------------------------------------------------------
        % Options:

        % ------------------------------------------------------------------------------
        % The standard methods
        % ------------------------------------------------------------------------------
        % cfg.removeNoise = '6P'
            % Six motion parameters from realignment
            % totaling 6 noise regressors
        % cfg.removeNoise = '6P+2P'
            % Six motion parameters from realignment
            % as well as mean WM/CSF
            % totaling 8 noise regressors
        % cfg.removeNoise = '6P+2P+GSR'
            % Six motion parameters from realignment
            % as well as GSR
            % totaling 9 noise regressors

        % ------------------------------------------------------------------------------
        % The expanded methods
        % ------------------------------------------------------------------------------
        % cfg.removeNoise = '24P'
            % Involves regressing out the following:
            % 1) 6 motion parameters from realignment
            % 3) the temporal derivatives of the above (calculated as the backwards difference)
            % 4) the squares of the above
            % totaling 24 noise regressors. (6*2*2) 
        % cfg.removeNoise = '24P+8P'
            % Involves regressing out the following:
            % 1) 6 motion parameters from realignment
            % 2) mean WM and CSF signals
            % 3) the temporal derivatives of the above (calculated as the backwards difference)
            % 4) the squares of the above
            % totaling 32 noise regressors. (6*2*2) + (2*2*2)
        % cfg.removeNoise = '24P+8P+4GSR'
            % 24P+8P as well as 4GSR (expanded)
            % totaling 36 noise regressors. (6*2*2) + (2*2*2) + (1*2*2)           
        % cfg.removeNoise = '24P+8P+SpikeReg'
            % 24P+8P as well as spike regression calculated by thresholding the output from GetFDJenk.m (as per Ciric/Satterthwaite)
            % totaling 24P+8P + N noise regressors (where N varies per subject depending on the number of suprathreshold fd spikes)
        % cfg.removeNoise = '24P+8P+4GSR+SpikeReg'
            % The 24P+8P+4GSR model with spike regression
            % AKA, the 24P+8P+SpikeReg model with 4GSR
            % So, 36 params + N spike regressors
            % Same as in Satterthwaite

        % Extra, not in paper
        % cfg.removeNoise = '24P+4GSR'
            % 24P as well as 4GSR (expanded)
            % totaling 28 noise regressors. (6*2*2) + (1*2*2)
        % ------------------------------------------------------------------------------
            
        % ------------------------------------------------------------------------------
        % The aCompCor methods.
        % ------------------------------------------------------------------------------
        % cfg.removeNoise = '12P+aCC'
            % As per Muschelli et al 2014
            % 1) 6 motion parameters from realignment
            % 2) the temporal derivatives of the above (calculated as the backwards difference)
            % 3) 10 regressors from aCompCor
                % 5 from WM and 5 from CSF
            % totaling 22 regressors. (6*2) + 5 + 5
        % cfg.removeNoise = '12P+aCC'
            % As per Muschelli et al 2014 but using enough PCs to account for 50% of the variance
        % cfg.removeNoise = '24P+aCC'
            % 1) 6 motion parameters from realignment
            % 2) the temporal derivatives of the above (calculated as the backwards difference)
            % 3) the squares of the above
            % 4) 10 regressors from aCompCor
                % 5 from WM and 5 from CSF
            % totaling 34 regressors. (6*2*2) + 5 + 5
        % cfg.removeNoise = '24P+aCC+4GSR'
            % 24P+aCC as well as 4GSR (including temporal derivatives and squares)
            % totaling 38 noise regressors. (6*2*2) + 5 + 5 + 4        
        % cfg.removeNoise = '24P+aCC+SpikeReg'
            % 24P+aCC as well as spike regression calculated by thresholding the output from GetFDJenk.m (as per Ciric/Satterthwaite)
            % totaling 24P+aCC + N noise regressors (where N varies per subject depending on the number of suprathreshold fd spikes)
        % cfg.removeNoise = '24P+aCC+4GSR+SpikeReg'
            % 24P+aCC+4GSR and spike regression
            % totaling 38 + N noise regressors. (6*2*2) + 5 + 5 + 4 + N spikes
        % ------------------------------------------------------------------------------

        % ------------------------------------------------------------------------------
        % The ICA methods
        % ------------------------------------------------------------------------------
        % cfg.removeNoise = 'ICA-AROMA+2P'
            % ICA-AROMA as well as mean WM/CSF
        % cfg.removeNoise = 'ICA-AROMA+GSR'
            % ICA-AROMA as well as GSR
        % cfg.removeNoise = 'ICA-AROMA+2P+GSR'
            % ICA-AROMA as well as mean WM/CSF and GSR

        % cfg.removeNoise = 'ICA-AROMA+8P'
            % ICA-AROMA as well as mean expanded WM/CSF
        % cfg.removeNoise = 'ICA-AROMA+4GSR'
            % ICA-AROMA as well as 4GSR
        % cfg.removeNoise = 'ICA-AROMA+8P+4GSR'
            % ICA-AROMA as well as expanded WM/CSF and 4GSR

    fprintf('\n\t\t ----- Noise correction ----- \n\n');

    fprintf(1, '\t\t Subject: %s \n', cfg.subject)
    fprintf(1, '\t\t Noise correction: %s%s \n', cfg.removeNoise, cfg.suffix);
    fprintf(1, '\t\t Input file: %s \n', cfg.CleanIn);

    % ------------------------------------------------------------------------------
    % Setup output directory
    % ------------------------------------------------------------------------------
        % If the first prefix on the cfg.CleanIn file is an 's' 
        % then assume that the file has been smoothed during prepro_base
        % and set outdir to include the s prefix to denote running noise correction 
        % on smoothed data.
        if cfg.CleanIn(1) == 's'
            outdir = [cfg.preprodir,'s',cfg.removeNoise,'/']; 
        else
            outdir = [cfg.preprodir,cfg.removeNoise,'/']; 
        end

        if ~isempty(cfg.suffix)
            outdir = [outdir(1:end-1),cfg.suffix,'/'];
        end

        if exist(outdir) == 0
            fprintf(1,'\n\t\t Initialising outdir\n')
            mkdir(outdir)
        elseif exist(outdir) == 7
            fprintf(1,'\n\t\t Cleaning and re-initialising outdir\n')
            rmdir(outdir,'s')
            mkdir(outdir)
        end
       
        cd(outdir)

    % ------------------------------------------------------------------------------
    % Decompress input files if necessary
    % Note, due to storage constraints on MASSIVE, all output files are compressed
    % after completing on run_prepro.m.
    % This code was added in cases where I wanted to run additional pipelines after
    % the initial run_prepro job.
    % Since this prepro_noise.m expects .nii extension files (as per outputs of prepro_base.m),
    % I added this code to check whether the input files are compressed (.gz) or not.
    % If they are, decompress.
    % They will all be recompressed once again after run_prepro.m completes
    % ------------------------------------------------------------------------------
        filesToCheck = {cfg.CleanIn, ...
                        cfg.NuisanceIn_wm, ...
                        cfg.NuisanceIn_csf, ...
                        cfg.wmmask, ...
                        cfg.csfmask, ...
                        cfg.BrainMask};

        dirsOfFiles = {cfg.preprodir, ...
                        cfg.preprodir, ...
                        cfg.preprodir, ...
                        cfg.t1dir, ...
                        cfg.t1dir, ...
                        cfg.preprodir};

        for i = 1:length(filesToCheck)
            if exist([dirsOfFiles{i},filesToCheck{i},'.gz']) == 2
                fprintf(1, '\t\t Decompressing %s\n',filesToCheck{i});
                gunzip([dirsOfFiles{i},filesToCheck{i},'.gz'],dirsOfFiles{i})
                delete([dirsOfFiles{i},filesToCheck{i},'.gz'])
            elseif exist([dirsOfFiles{i},filesToCheck{i},'.gz']) == 0
                fprintf(1, '\t\t No need to decompress %s\n',filesToCheck{i});
            end
        end


    % ------------------------------------------------------------------------------
    % Setup switches
    % ------------------------------------------------------------------------------
        cfg.removeNoiseSplit = strsplit(cfg.removeNoise,'+');

        % ICA-AROMA
        if any(strmatch('ICA-AROMA',cfg.removeNoiseSplit,'exact')) == 1
            runICA = 1;
        else
            runICA = 0;
        end

        % Physiological noise
        if any(strmatch('2P',cfg.removeNoiseSplit,'exact')) == 1 || ...
            any(strmatch('8P',cfg.removeNoiseSplit,'exact')) == 1
            runPhys = 1;
        else
            runPhys = 0;
        end

        % GSR
        if any(strmatch('GSR',cfg.removeNoiseSplit,'exact')) == 1 || ...
            any(strmatch('4GSR',cfg.removeNoiseSplit,'exact')) == 1
            runGSR = 1;
        else
            runGSR = 0;
        end

        % aCompCor
        if any(strmatch('aCC',cfg.removeNoiseSplit,'exact')) == 1
            runaCC = 1;
            aCCversion = 5;
        elseif any(strmatch('aCC50',cfg.removeNoiseSplit,'exact')) == 1
            runaCC = 1;
            aCCversion = 50;
        else
            runaCC = 0;
        end

        % Spike regression
        if any(strmatch('SpikeReg',cfg.removeNoiseSplit,'exact')) == 1
            runSpikeReg = 1;
        else
            runSpikeReg = 0;
        end

    % ------------------------------------------------------------------------------
    % 0) ICA AROMA
    % ------------------------------------------------------------------------------
        if runICA == 1
            % set FSL output to nifti_gz because ICA-AROMA requires it!!!!
            setenv('FSLOUTPUTTYPE','NIFTI_GZ');

            scriptdir_ICA = '/gpfs/M2Home/projects/Monash076/Linden/scripts/Software/ICA-AROMA-master/';

            str = ['python2.7 ',scriptdir_ICA,'ICA_AROMA.py -in ',cfg.preprodir,cfg.CleanIn,' -out ',outdir,' -mc ',cfg.preprodir,'rp*.txt -m ',cfg.preprodir,cfg.BrainMask,' -tr ',num2str(cfg.TR)];
            system(str);
    
            % unzip
            system('gunzip -rf denoised_func_data_nonaggr.nii.gz');

            fprintf(1, '\n\t\t !!!! Overriding inputs for ICA-AROMA !!!! \n\n');
            % redefine clean in
            CleanInNew = ['ica_',cfg.CleanIn];
            cfg.CleanIn = CleanInNew;
            % redefine wm/csf files
            cfg.NuisanceIn_wm = CleanInNew;
            cfg.NuisanceIn_csf = CleanInNew;

            % Move ICA-AROMA output file up a step in the directory and rename
            movefile('denoised_func_data_nonaggr.nii',[cfg.preprodir,cfg.CleanIn])

            % Set FSL output back to nifti
            setenv('FSLOUTPUTTYPE','NIFTI');
        end

    % ------------------------------------------------------------------------------
    % 1) Motion parameters
    % ------------------------------------------------------------------------------
        % read in motion
        mfile = dir([cfg.preprodir,'rp*txt']);
        mov = dlmread([cfg.preprodir,mfile(1).name]);

    % ------------------------------------------------------------------------------
    % 2) Physiological time series
    % ------------------------------------------------------------------------------
        if runPhys == 1
            % Extract time series

            % White Matter
            str = [cfg.fsldir,'fslmeants -i ',cfg.preprodir,cfg.NuisanceIn_wm,' -o wmTS.txt -m ',cfg.t1dir,cfg.wmmask];
            system(str); 
            wmTS = dlmread('wmTS.txt');

            % CSF
            str = [cfg.fsldir,'fslmeants -i ',cfg.preprodir,cfg.NuisanceIn_csf,' -o csfTS.txt -m ',cfg.t1dir,cfg.csfmask];
            system(str); 
            csfTS = dlmread('csfTS.txt');
        end

    % ------------------------------------------------------------------------------
    % 2.1) Global Signal Regression
    % ------------------------------------------------------------------------------
        if runGSR == 1
            % Global
            str = [cfg.fsldir,'fslmeants -i ',cfg.preprodir,cfg.CleanIn,' -o gsTS.txt -m ',cfg.preprodir,cfg.BrainMask];
            system(str); 
            gsTS = dlmread('gsTS.txt');
        end

    % ------------------------------------------------------------------------------
    % 3) CompCor
    % ------------------------------------------------------------------------------
        if runaCC == 1
            % 1) White Matter
            if aCCversion == 5
                fprintf(1, 'Running aCompCor: white matter. 5 PCs. \n');
            elseif aCCversion == 50
                fprintf(1, 'Running aCompCor: white matter. 50%% PCs. \n');
            end

            % Check if tissue map is empty
            [~,temp] = read([cfg.t1dir,cfg.wmmask]);
            if numel(unique(temp)) == 1
                fprintf(1, 'WARNING! Your tissue mask is empty. Please check.\n');
                return
            else
                % Extract nuisance time course
                str = [cfg.fsldir,'fslmeants -i ',cfg.preprodir,cfg.NuisanceIn_wm,' -o wmTS_aCC.txt -m ',cfg.t1dir,cfg.wmmask,' --showall'];
                system(str); 

                % Read in cfg.wm time courses
                wmTS = dlmread('wmTS_aCC.txt');
                wmTS(1:3,:) = []; % delete first 3 rows, which list voxel coords

                % run cfg.wm aCompCor
                if aCCversion == 5
                    wm_aCompCor = CompCor(wmTS,'retain',[],5);
                elseif aCCversion == 50
                    wm_aCompCor = CompCor(wmTS,'50pc');
                end

                clear wmTS % clear to save memory
            end
            clear temp

            % 2) CSF
            if aCCversion == 5
                fprintf(1, 'Running aCompCor: csf. 5 PCs. \n');
            elseif aCCversion == 50
                fprintf(1, 'Running aCompCor: csf. 50%% PCs. \n');
            end

            % Check if tissue map is empty
            [~,temp] = read([cfg.t1dir,cfg.csfmask]);
            if numel(unique(temp)) == 1
                fprintf(1, 'WARNING! Your tissue mask is empty. Please check.\n');
                return
            else
                % Extract nuisance time course
                str = [cfg.fsldir,'fslmeants -i ',cfg.preprodir,cfg.NuisanceIn_csf,' -o csfTS_aCC.txt -m ',cfg.t1dir,cfg.csfmask,' --showall'];
                system(str); 

                % Read in cfg.csf time courses
                csfTS = dlmread('csfTS_aCC.txt');
                csfTS(1:3,:) = []; % delete first 3 rows, which list voxel coords

                % run cfg.csf aCompCor
                if aCCversion == 5
                    csf_aCompCor = CompCor(csfTS,'retain',[],5);
                elseif aCCversion == 50
                    csf_aCompCor = CompCor(csfTS,'50pc');
                end

                clear csfTS % clear to save memory
            end
            clear temp
        end
    
    % ------------------------------------------------------------------------------
    % 4) Spike regression
    % ------------------------------------------------------------------------------
        if runSpikeReg == 1
            % calculate Jenkinson FD
            fd = GetFDJenk(mov);

            % generate spike regressors
            spikereg = GetSpikeRegressors(fd,0.25);

            clear fd
        end
    
    % ------------------------------------------------------------------------------
    % Generate noiseTS variable for fsl_regfilt
    % ------------------------------------------------------------------------------
        fprintf(1,'\n\t\t ----- Generating noiseTS ----- \n\n');

        % If ICA wasn't run, then start noiseTS with movement params from realignment
        if runICA == 0
            fprintf(1,'\t\t Adding motion params \n');
            noiseTS = mov;
        % If it was run, then we dont regress out mov params
        elseif runICA == 1
            fprintf(1,'\t\t Skipping motion params \n');
            noiseTS = [];
        end

        % concatenate phys time series
        if runPhys == 1
            fprintf(1,'\t\t Adding phys signals \n');
            noiseTS = [noiseTS wmTS csfTS];
        end

        % concatenate global time series
        if runGSR == 1
            fprintf(1,'\t\t Adding global signal \n');
            noiseTS = [noiseTS gsTS];
        end

        % If expanded model is selected, get derivatives and squares
        if any(strmatch('24P',cfg.removeNoiseSplit,'exact')) == 1 | ...
            any(strmatch('8P',cfg.removeNoiseSplit,'exact')) == 1 | ...
            any(strmatch('4GSR',cfg.removeNoiseSplit,'exact')) == 1
            fprintf(1,'\t\t Getting derivatives and squares \n');
            noiseTS = GetDerivatives(noiseTS);
        elseif any(strmatch('12P',cfg.removeNoiseSplit,'exact')) == 1
            fprintf(1,'\t\t Getting derivatives only \n');
            noiseTS = GetDerivatives(noiseTS,0);
        else
            % Otherwise, just detrend (note, detrending done as part og GetDerivatives.m)
            fprintf(1,'\t\t Detrending only \n');
            noiseTS = detrend(noiseTS,'linear');
        end

        % We dont get expansion terms on spike regressors or aCompCor
        % concatenate the spike regressors and compcor stuff after derivatives.
        
        % concatenate aCompCor
        if runaCC == 1
            fprintf(1,'\t\t Adding aCompCor signals \n');
            noiseTS = [noiseTS wm_aCompCor csf_aCompCor];
        end

        % concatenate spike regressors
        if runSpikeReg == 1
            fprintf(1,'\t\t Adding spike regressors \n');
            noiseTS = [noiseTS spikereg];
        end

    % ------------------------------------------------------------------------------
    % Clean data: fsl_regfilt
    % ------------------------------------------------------------------------------
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
        if runICA == 1
            % Move back
            movefile([cfg.preprodir,cfg.CleanIn],outdir)
        end

    % ------------------------------------------------------------------------------
    % Bandpass filter with REST
    % ------------------------------------------------------------------------------
        fprintf(1,'\n\t\t ----- Running bandpass filtering ----- \n\n');

        FiltIn = 'epi_clean.nii';
        CUTNUMBER = 1;

        cd(outdir)
        % creater input directory for REST function
        mkdir('temp');
        filtdir = [outdir,'temp/'];
        movefile(FiltIn,filtdir,'f');
        
        switch cfg.WhichNii
            case '4D'
                % bandpass epis using rest
                rest_bandpass(filtdir,cfg.TR,cfg.LowPass,cfg.HighPass,'No',[cfg.preprodir,cfg.BrainMask],CUTNUMBER)

                % Move output files back to outdir
                movefile('temp/*',outdir)
                movefile('temp_filtered/*',outdir)
                rmdir('temp','s')
                rmdir('temp_filtered','s')

                % Rename
                movefile('Filtered_4DVolume.nii','epi_prepro.nii')

            case '3D'
                % split EPI into 3D files
                cd(filtdir)
                spm_file_split(FiltIn)

                % Move 4D file back to outdir
                movefile(FiltIn,outdir,'f')

                cd(outdir)
                % bandpass epis using rest
                rest_bandpass(filtdir,cfg.TR,cfg.LowPass,cfg.HighPass,'No',[cfg.preprodir,cfg.BrainMask],CUTNUMBER)

                % Move output files back to outdir
                rmdir(filtdir,'s')
                movefile('temp_filtered/*',outdir)
                rmdir('temp_filtered','s')

                % Rename
                movefile('Filtered_4DVolume.nii','epi_prepro.nii')
        end

    % ------------------------------------------------------------------------------
    % Spatially smooth the data
    % ------------------------------------------------------------------------------
        switch cfg.smoothing
            case 'after'
                cd(outdir)

                fprintf('\n\t\t ----- Spatial smoothing ----- \n\n');

                SmoothIn = 'epi_prepro.nii';

                str = [cfg.afnidir,'3dBlurInMask -input ',SmoothIn,' -FWHM ',num2str(cfg.kernel),' -mask ',cfg.preprodir,cfg.BrainMask,' -prefix smoothed'];
                system(str);
                % convert to nifti
                system([cfg.afnidir,'3dAFNItoNIFTI smoothed+tlrc']);
                % delete afni outputs
                delete('smoothed+tlrc*')
                % rename output file
                movefile('smoothed.nii',['s',SmoothIn])
            case 'none'
                fprintf('\n\t\t !!!! Skipping spatial smoothing !!!! \n\n');
                % do nothing
        end

    fprintf('\n\t\t ----- Noise correction complete ----- \n\n');
end
