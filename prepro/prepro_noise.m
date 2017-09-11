function [noiseTS,outdir,noiseTSz] = prepro_noise(cfg)
    % 1 - Correct noise in EPI output from step 11 of prepro_base script (or output from step 12 in case of ICA-AROMA)
    %       See below for details
    % 2 - Bandpass filter (includes demeaning)
    % 3 - Spatial smoothing (not in case of ICA-AROMA)
    % 
    % Copyright (C) 2017, Linden Parkes <lindenparkes@gmail.com>,
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
    fprintf(1, '\t\t Noise correction: %s \n', cfg.removeNoise);
    fprintf(1, '\t\t Input file: %s \n', cfg.CleanIn);

    if cfg.runReg == 1
        fprintf(1, '\t\t Nuisance regression: extract and regress \n');
    elseif cfg.runReg == 0
        fprintf(1, '\t\t Nuisance regression: extract only \n');
    end

    fprintf('\n\t\t ---------------------------- \n\n');

    % ------------------------------------------------------------------------------
    % Setup output directory
    % ------------------------------------------------------------------------------
        outdir = [cfg.preprodir,cfg.removeNoise,'/']; 

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
            any(strmatch('4P',cfg.removeNoiseSplit,'exact')) == 1 || ...
            any(strmatch('8P',cfg.removeNoiseSplit,'exact')) == 1
            runPhys = 1;
        else
            runPhys = 0;
        end

        % GSR
        if any(strmatch('GSR',cfg.removeNoiseSplit,'exact')) == 1 || ...
            any(strmatch('2GSR',cfg.removeNoiseSplit,'exact')) == 1 || ...
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

        % Power 2014 scrubbing
        if any(strmatch('JP14Scrub',cfg.removeNoiseSplit,'exact')) == 1
            runJP14Scrub = 1;
        else
            runJP14Scrub = 0;
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

            if runJP14Scrub == 1
                ICA_outdir = [cfg.preprodir,'JP14_ICA-AROMA_output/'];
            elseif runJP14Scrub == 0
                ICA_outdir = [cfg.preprodir,'ICA-AROMA_output/'];
            end

            % First check if ICA-AROMA has been run already
            if exist(ICA_outdir) == 0
                % If it doesn't, run ICA-AROMA
                
                % set FSL output to nifti_gz because ICA-AROMA requires it!!!!
                setenv('FSLOUTPUTTYPE','NIFTI_GZ');

                fprintf(1,'\n\t\t Initialising ICA_outdir\n')
                mkdir(ICA_outdir)
                cd(ICA_outdir)
    
                if cfg.runReg == 1
                    str = ['python2.7 ',cfg.scriptdir_ICA,'ICA_AROMA.py -in ',cfg.preprodir,cfg.CleanIn,' -out ',ICA_outdir,' -mc ',cfg.preprodir,'raw_mov/rp*.txt -m ',cfg.preprodir,cfg.BrainMask,' -tr ',num2str(cfg.TR)];
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

                    % rename output file
                    movefile('denoised_func_data_nonaggr.nii',[cfg.preprodir,cfg.CleanIn])

                elseif cfg.runReg == 0
                    str = ['python2.7 ',cfg.scriptdir_ICA,'ICA_AROMA.py -in ',cfg.preprodir,cfg.CleanIn,' -out ',ICA_outdir,' -mc ',cfg.preprodir,'raw_mov/rp*.txt -m ',cfg.preprodir,cfg.BrainMask,' -tr ',num2str(cfg.TR),' -den no'];
                    system(str);
                end

                % Set FSL output back to nifti
                setenv('FSLOUTPUTTYPE','NIFTI');
            elseif exist(ICA_outdir) == 7
                fprintf(1,'\n\t\t ICA-AROMA has already been done: skipping \n')
                cd(ICA_outdir)

                fprintf(1, '\n\t\t !!!! Overriding inputs for ICA-AROMA !!!! \n\n');
                % redefine clean in
                CleanInNew = ['ica_',cfg.CleanIn];
                cfg.CleanIn = CleanInNew;
                % redefine wm/csf files
                cfg.NuisanceIn_wm = CleanInNew;
                cfg.NuisanceIn_csf = CleanInNew;

                % move ICA cleaned file
                movefile(cfg.CleanIn,cfg.preprodir)
            end
            % Change back to primary outdir
            cd(outdir)
        end

    % ------------------------------------------------------------------------------
    % 1) Motion parameters
    % ------------------------------------------------------------------------------
        % read in motion
        mfile = dir([cfg.preprodir,'raw_mov/rp*.txt']);
        mov = dlmread([cfg.preprodir,'raw_mov/',mfile(1).name]);

    % ------------------------------------------------------------------------------
    % 2) Physiological time series
    % ------------------------------------------------------------------------------
        if runPhys == 1
            % Extract time series

            % White Matter
            str = [cfg.fsldir,'fslmeants -i ',cfg.preprodir,cfg.NuisanceIn_wm,' -o wmTS.txt -m ',cfg.t1dir,cfg.wmmask];
            system(str); 
            wmTS = dlmread('wmTS.txt');

            % also get the other less eroded wm masks.
            % Note, this code is a bit lazy and will end up with a duplicate .txt file.
            % e.g., wmTS.txt is probablity be identical to wm_e5_TS.txt
            for i = 1:length(cfg.wm)
                [hdr,data] = read([cfg.t1dir,cfg.wm{i}]);
                if any(data(:)) == 1;
                    str = [cfg.fsldir,'fslmeants -i ',cfg.preprodir,cfg.NuisanceIn_wm,' -o wm_e',num2str(i),'_TS.txt -m ',cfg.t1dir,cfg.wm{i}];
                    system(str); 
                end
            end                

            % CSF
            str = [cfg.fsldir,'fslmeants -i ',cfg.preprodir,cfg.NuisanceIn_csf,' -o csfTS.txt -m ',cfg.t1dir,cfg.csfmask];
            system(str); 
            csfTS = dlmread('csfTS.txt');

            for i = 1:length(cfg.csf)
                [hdr,data] = read([cfg.t1dir,cfg.csf{i}]);
                if any(data(:)) == 1;
                    str = [cfg.fsldir,'fslmeants -i ',cfg.preprodir,cfg.NuisanceIn_csf,' -o csf_e',num2str(i),'_TS.txt -m ',cfg.t1dir,cfg.csf{i}];
                    system(str); 
                end
            end

            % Gray Matter
            % Note, not used for denoising.
            % We only retain a pre-cleaned GM signal so we can examinine wm-gm correlations as a function of wm erosion (Power et al., 2017)
            % Note correlation to WM will likely be 0.1-0.2 higher than in Power et al., because Power uses the ribbon only to get GM (ribbon is better)
            str = [cfg.fsldir,'fslmeants -i ',cfg.preprodir,cfg.CleanIn,' -o gmTS.txt -m ',cfg.t1dir,cfg.gmmask];
            system(str);             
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
            fprintf(1,'\n\t\t ----- Anatomical CompCor ----- \n\n')

            % 1) White Matter
            if aCCversion == 5
                fprintf(1, '\t\t Running aCompCor: white matter. 5 PCs. \n');
            elseif aCCversion == 50
                fprintf(1, '\t\t Running aCompCor: white matter. 50%% PCs. \n');
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

                % Save out number of components
                aCC_num_wm = size(wm_aCompCor,2);
                dlmwrite('aCC_num_wm.txt',aCC_num_wm)

                clear wmTS % clear to save memory
            end
            clear temp

            % 2) CSF
            if aCCversion == 5
                fprintf(1, '\t\t Running aCompCor: csf. 5 PCs. \n');
            elseif aCCversion == 50
                fprintf(1, '\t\t Running aCompCor: csf. 50%% PCs. \n');
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

                % Save out number of components
                aCC_num_csf = size(csf_aCompCor,2);
                dlmwrite('aCC_num_csf.txt',aCC_num_csf)

                clear csfTS % clear to save memory
            end
            clear temp
        end
    
    % ------------------------------------------------------------------------------
    % 4) Spike regression
    % ------------------------------------------------------------------------------
        if runSpikeReg == 1
            fprintf(1,'\n\t\t ----- Generating spike regressors ----- \n\n');

            % calculate Jenkinson FD
            fd = GetFDJenk(mov);

            % generate spike regressors
            spikereg = GetSpikeRegressors(fd,0.25);

            clear fd
        end
    
    % ------------------------------------------------------------------------------
    % Generate noiseTS
    % ------------------------------------------------------------------------------
        fprintf(1,'\n\t\t ----- Generating noiseTS ----- \n\n');

        % If ICA wasn't run, then start noiseTS with movement params from realignment
        if runICA == 0
            fprintf(1,'\t\t Adding motion params \n');

            % Get expansions for motion
            if any(strmatch('24P',cfg.removeNoiseSplit,'exact')) == 1
                fprintf(1,'\t\t Getting derivatives and squares \n');
                noiseTS = GetDerivatives(mov);
            elseif any(strmatch('12P',cfg.removeNoiseSplit,'exact')) == 1
                fprintf(1,'\t\t Getting derivatives only \n');
                noiseTS = GetDerivatives(mov,0);
            else
                % Otherwise, just detrend (note, detrending done as part of GetDerivatives.m)
                fprintf(1,'\t\t Detrending only \n');
                noiseTS = detrend(mov);                
            end
        % If it was run, then we dont regress out mov params
        elseif runICA == 1
            fprintf(1,'\t\t Skipping motion params \n');
            noiseTS = [];
        end

        % concatenate phys time series
        if runPhys == 1
            fprintf(1,'\t\t Adding phys signals \n');
            physTS = [wmTS csfTS];

            if any(strmatch('8P',cfg.removeNoiseSplit,'exact')) == 1
                fprintf(1,'\t\t Getting derivatives and squares \n');
                physTS = GetDerivatives(physTS);
            elseif any(strmatch('4P',cfg.removeNoiseSplit,'exact')) == 1
                fprintf(1,'\t\t Getting derivatives only \n');
                physTS = GetDerivatives(physTS,0);
            else
                % Otherwise, just detrend (note, detrending done as part of GetDerivatives.m)
                fprintf(1,'\t\t Detrending only \n');
                physTS = detrend(physTS);                
            end

            noiseTS = [noiseTS physTS];
        end

        % concatenate global time series
        if runGSR == 1
            fprintf(1,'\t\t Adding global signal \n');

            if any(strmatch('4GSR',cfg.removeNoiseSplit,'exact')) == 1
                fprintf(1,'\t\t Getting derivatives and squares \n');
                gsTS = GetDerivatives(gsTS);
            elseif any(strmatch('2GSR',cfg.removeNoiseSplit,'exact')) == 1
                fprintf(1,'\t\t Getting derivatives only \n');
                gsTS = GetDerivatives(gsTS,0);
            else
                % Otherwise, just detrend (note, detrending done as part of GetDerivatives.m)
                fprintf(1,'\t\t Detrending only \n');
                gsTS = detrend(gsTS);                
            end

            noiseTS = [noiseTS gsTS];
        end

        % We dont get expansion terms on spike regressors or aCompCor
        
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
    % Nuisance regression
    % ------------------------------------------------------------------------------
        if cfg.runReg == 1

            [hdr,data] = read([cfg.preprodir,cfg.CleanIn]);

            % reshape to 2d
            dim = size(data);
            data = reshape(data,[],dim(4));

            if runJP14Scrub == 1
                fprintf(1,'\n\t\t ----- Running nuisance regression with scrubbing ----- \n\n');
                % Read in scrubbing mask
                scrubmask = logical(dlmread([cfg.preprodir,'JP14_ScrubMask.txt']));
                [data_out zb noiseTSz] = JP14_regress_nuisance(data,noiseTS,~scrubmask(:,2));
            elseif runJP14Scrub == 0
                fprintf(1,'\n\t\t ----- Running nuisance regression without scrubbing ----- \n\n');
                [data_out zb noiseTSz] = JP14_regress_nuisance(data,noiseTS);
            end

            CleanOut = 'epi_clean.nii';
            data_out = reshape(data_out,dim);
            write(hdr,data_out,CleanOut)            
            
            % write out noiseTS incase people want to model nuisance at SPM
            dlmwrite('noiseTS.txt',noiseTS,'delimiter','\t','precision','%.6f');
            dlmwrite('noiseTSz.txt',noiseTSz,'delimiter','\t','precision','%.6f');
        elseif cfg.runReg == 0
            fprintf(1,'\n\t\t ----- Skipping nuisance regression ----- \n\n');
        end

    % ------------------------------------------------------------------------------
    % Bandpass filter with REST
    % ------------------------------------------------------------------------------
        if cfg.runReg == 1
            FiltIn = CleanOut;
        elseif cfg.runReg == 0
            FiltIn = cfg.CleanIn;
            copyfile([cfg.preprodir,FiltIn],outdir)
        end

        % load nifti
        [hdr,data] = read(FiltIn);

        % reshape to 2d
        dim = size(data);
        data = reshape(data,[],dim(4));
        % Put time on first dimension
        data = data';
        numVoxels = size(data,2);

        % load in brain mask
        [hdr_mask,data_mask] = read([cfg.preprodir,cfg.BrainMask]);
        data_mask = reshape(data_mask,[],1);
        data_mask = logical(data_mask);

        % mask out non-brain voxels
        % Note, we only do this to speed up the interpolation step that is run only for JP14Scrub method
        data = data(:,data_mask);

        if runJP14Scrub == 1
            fprintf(1,'\n\t\t ----- Running bandpass filtering with interpolation ----- \n\n');
            % Read in scrubbing mask
            scrubmask = logical(dlmread([cfg.preprodir,'JP14_ScrubMask.txt']));
            TRtimes = ([1:cfg.tN]')*cfg.TR;

            voxbinsize = 500;
            voxbin = 1:voxbinsize:size(data,2);
            voxbin = [voxbin size(data,2)];

            data_surrogate = zeros(size(data,1),size(data,2));

            for v = 1:numel(voxbin)-1 % this takes huge RAM if all voxels
                fprintf(1, 'Processing voxel bin %u/%u\n', v,numel(voxbin)-1);
                data_surrogate(:,voxbin(v):voxbin(v+1)) = JP14_getTransform(data(:,voxbin(v):voxbin(v+1)),TRtimes,scrubmask(:,2));
            end

            % insert surrogate data into real data at censored time points
            data(scrubmask(:,2),:) = data_surrogate(scrubmask(:,2),:);
        elseif runJP14Scrub == 0
            fprintf(1,'\n\t\t ----- Running bandpass filtering without interpolation ----- \n\n');
        end

        % bandpass filter the masked data
        data_filtered = rest_IdealFilter(data, cfg.TR, [cfg.HighPass, cfg.LowPass]);

        % put filtered data back with non-brain voxels 
        data_out = zeros(cfg.tN,numVoxels);
        data_out(:,data_mask) = data_filtered;

        % reshape back to 4d
        data_out = data_out';
        data_out = reshape(data_out,dim);

        % write nifti
        FiltOut = 'epi_filtered.nii';
        write(hdr,data_out,FiltOut)

    % ------------------------------------------------------------------------------
    % Spatially smooth the data
    % ------------------------------------------------------------------------------
        switch cfg.smoothing
            case 'after'
                fprintf('\n\t\t ----- Spatial smoothing ----- \n\n');

                SmoothIn = FiltOut;
                SmoothEPI(SmoothIn,cfg.kernel,cfg.tN)
                movefile(['s',SmoothIn],'epi_prepro.nii')

                % str = [cfg.afnidir,'3dBlurInMask -input ',SmoothIn,' -FWHM ',num2str(cfg.kernel),' -mask ',cfg.preprodir,cfg.BrainMask,' -prefix smoothed'];
                % system(str);
                % % convert to nifti
                % system([cfg.afnidir,'3dAFNItoNIFTI smoothed*']);
                % % rename output file
                % movefile('smoothed.nii','epi_prepro.nii')

                % % delete afni outputs
                % delete('*.BRIK')
                % delete('*.HEAD')

            case 'none'
                fprintf('\n\t\t !!!! Skipping spatial smoothing !!!! \n\n');
                % Simply rename the filtered epi to epi_prepro.nii so that the final output file is always the same file name
                movefile(FiltOut,'epi_prepro.nii')
            case 'before'
                movefile(FiltOut,'epi_prepro.nii')
        end

    % ------------------------------------------------------------------------------
    % Put ICA-AROMA output back in ICA dir
    % ------------------------------------------------------------------------------
        if runICA == 1
            % Move back
            movefile([cfg.preprodir,cfg.CleanIn],ICA_outdir)
        end

    fprintf('\n\t\t ----- Noise correction complete ----- \n\n');
end
