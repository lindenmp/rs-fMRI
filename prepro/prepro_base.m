function [tN,gm,wm,csf,epiBrainMask,t1BrainMask,BrainMask,gmmask,wmmask,csfmask,dvars,dvarsExtract,fdThr,dvarsThr,exclude,outEPI] = prepro_base(cfg)
    % This script performs the base preprocessing, upon all variants of nuissance removal is run
    % These steps includes:
    %
    % ------------------------------------------------------------------------------
    % T1 stream
    % ------------------------------------------------------------------------------
    % 1 - Segment T1 using SPM (Generate tissue masks)

    % ------------------------------------------------------------------------------
    % EPI stream
    % ------------------------------------------------------------------------------
    % 2 - Discard first 4 volumes (optional. set using input argument)
    % 3 - Slice-timing correction (optional. set using input argument)
    % 4 - Despiking (optional. set using input argument)
    % 5 - EPI Realignment
    % 6 - Co-registration of realigned images to native T1
    % 7 - Spatially normalize T1 to MNI template
    % 8 - Application of T1-spatial normalization parameters to coregistered EPI and tissue masks
    % 9 - Mask out non-brain voxels
    % 10 - Linear detrending of realigned EPI time series (uses REST) (optional. set using input argument)
    % 11 - Modal intensity normalisation (to 1000) (optional. set using input argument)
    % 12 - Spatial smoothing
    %
    % Copyright (C) 2017, Linden Parkes <lindenparkes@gmail.com>,
    %
    % ------------------------------------------------------------------------------

    fprintf('\n\t\t ----- Base preprocessing ----- \n\n');

    fprintf(1, '\t\t Subject: %s \n', cfg.subject)

    if cfg.discard == 1
        fprintf(1, '\t\t Discard EPI volumes: yes \n');
    elseif cfg.discard == 0
        fprintf(1, '\t\t Discard EPI volumes: no \n');
    end

    if cfg.slicetime == 1
        fprintf(1, '\t\t Slicetime EPI: yes \n');
    elseif cfg.slicetime == 0
        fprintf(1, '\t\t Slicetime EPI: no \n');
    end

    if cfg.despike == 1
        fprintf(1, '\t\t Despike EPI: yes \n');
    elseif cfg.despike == 0
        fprintf(1, '\t\t Despike EPI: no \n');
    end

    if cfg.intnorm == 1
        fprintf(1, '\t\t Normalise EPI intensity: yes \n');
    elseif cfg.intnorm == 0
        fprintf(1, '\t\t Normalise EPI intensity: no \n');
    end

    if cfg.detr == 1
        fprintf(1, '\t\t Detrend EPI: yes \n');
    elseif cfg.detr == 0
        fprintf(1, '\t\t Detrend EPI: no \n');
    end

    if cfg.meanback == 1
        fprintf(1, '\t\t\t Add mean back: yes \n');
    elseif cfg.meanback == 0
        fprintf(1, '\t\t\t Add mean back: no \n');
    end

    fprintf('\n\t\t ------------------------------ \n\n');

    % ------------------------------------------------------------------------------
    % Set prepro dir
    % ------------------------------------------------------------------------------
        if exist(cfg.preprodir) == 0
            fprintf(1,'\t\t Initialising preprodir\n')
            mkdir(cfg.preprodir)
        elseif exist(cfg.preprodir) == 7
            fprintf(1,'\t\t Cleaning and re-initialising preprodir\n')
            rmdir(cfg.preprodir,'s')
            mkdir(cfg.preprodir)
        end

    % ------------------------------------------------------------------------------
    % Decompress files
    % ------------------------------------------------------------------------------
        % t1
        cd(cfg.t1dir);
        % Check extension of cfg.t1name
        % If user entered in a file with .gz extension, remove it and only retain the .nii extension.
        % note, the next section of code will search for .gz file and decompress if necessary
        [~,name,ext] = fileparts(cfg.t1name);
        switch ext
            case '.gz'
                cfg.t1name = name;
        end
        
        if exist(cfg.t1name) == 0
            if exist([cfg.t1name,'.gz']) == 2
                runDecompt1 = 1;
            elseif exist([cfg.t1name,'.gz']) == 0
                fprintf(1, 'Warning: raw t1 not found!\n');
            end
        elseif exist(cfg.t1name) == 2;
            runDecompt1 = 0;
        end

        if runDecompt1 == 1
            fprintf(1, '\t\t Decompressing t1...\n');
            gunzip([cfg.t1name,'.gz'])
            delete([cfg.t1name,'.gz'])
        end

        % Extra t1 (conc_TMS_fMRI project only)
        if isfield(cfg, 't14norm')
            [~,name,ext] = fileparts(cfg.t14norm);
            switch ext
                case '.gz'
                    cfg.t14norm = name;
            end
        
            if exist(cfg.t14norm) == 0
                if exist([cfg.t14norm,'.gz']) == 2
                    runDecompt1 = 1;
                elseif exist([cfg.t14norm,'.gz']) == 0
                    fprintf(1, 'Warning: raw t14norm not found!\n');
                end
            elseif exist(cfg.t14norm) == 2;
                runDecompt1 = 0;
            end

            if runDecompt1 == 1
                fprintf(1, '\t\t Decompressing t14norm...\n');
                gunzip([cfg.t14norm,'.gz'])
                delete([cfg.t14norm,'.gz'])
            end
        end

        % EPI
        cd(cfg.rawdir);
        % Check extension of cfg.EPI
        % If user entered in a file with .gz extension, remove it and only retain the .nii extension.
        % note, the next section of code will search for .gz file and decompress if necessary
        [~,name,ext] = fileparts(cfg.EPI);
        switch ext
            case '.gz'
                cfg.EPI = name;
        end
        
        if exist(cfg.EPI) == 0
            if exist([cfg.EPI,'.gz']) == 2
                runDecompEPI = 1;
            elseif exist([cfg.EPI,'.gz']) == 0
                fprintf(1, 'Warning: raw EPI not found!\n');
            end
        elseif exist(cfg.EPI) == 2;
            runDecompEPI = 0;
        end

        if runDecompEPI == 1
            fprintf(1, '\t\t Decompressing epi...\n');
            gunzip([cfg.EPI,'.gz'])
            delete([cfg.EPI,'.gz'])
        end

    % ------------------------------------------------------------------------------
    % Preprocess T1
    % In this step, we take the subject's T1 image and segment it into the three
    % basic tissue types using SPM8 New Segment routine
    % The tissue types are:
    % 1) grey matter
    % 2) white matter
    % 3) cerebrospinal fluid
    % The results are probabilistic brain images where higher values in voxels
    % represent a greater chance that a given voxel is a given tissue type
    % ------------------------------------------------------------------------------
        fprintf(1, '\n\t\t ----- Processing T1 ----- \n\n');
        cd([cfg.datadir,cfg.subject])

        % Clean and reinitialise T1 dir
        movefile([cfg.t1dir,cfg.t1name],[cfg.datadir,cfg.subject])
        if ~isempty(dir([cfg.t1dir,'*.json'])) == 1
            fname = dir([cfg.t1dir,'*.json']);
            movefile([cfg.t1dir,fname.name],[cfg.datadir,cfg.subject])
        end

        if isfield(cfg, 't14norm')
            movefile([cfg.t1dir,cfg.t14norm,'*'],[cfg.datadir,cfg.subject])
        end

        delete([cfg.t1dir,'*'])

        movefile([cfg.datadir,cfg.subject,'/',cfg.t1name],cfg.t1dir)
        if ~isempty(dir([cfg.datadir,cfg.subject,'/*.json'])) == 1
            movefile([cfg.datadir,cfg.subject,'/',fname.name],cfg.t1dir)
        end
        clear fname

        if isfield(cfg, 't14norm')
            movefile([cfg.datadir,cfg.subject,'/',cfg.t14norm,'*'],cfg.t1dir)
        end
        
        cd(cfg.t1dir);

        % First crop out neck
        outname = ['c',cfg.t1name];
        system([cfg.fsldir,'robustfov -i ',cfg.t1name,' -r ',outname]);
        cfg.t1name = outname;

        % N4 Bias field correction w/ ANTs
        outname = ['N4',cfg.t1name];
        system([cfg.antsdir,'N4BiasFieldCorrection -d 3 -i ',cfg.t1name,' -c [100x100x100x100,0.000001] -o ',outname]);
        cfg.t1name = outname;

        % Do fDown T1 if it exists (conc_TMS_fMRI project only)
        if isfield(cfg, 't14norm')
            outname = ['c',cfg.t14norm];
            system([cfg.fsldir,'robustfov -i ',cfg.t14norm,' -r ',outname]);
            cfg.t14norm = outname;

            outname = ['N4',cfg.t14norm];
            system([cfg.antsdir,'N4BiasFieldCorrection -d 3 -i ',cfg.t14norm,' -c [100x100x100x100,0.000001] -o ',outname]);
            cfg.t14norm = outname;
        end

        % Tissue segment T1 with SPM
        SegmentT1([cfg.t1dir,cfg.t1name],cfg.spmdir,0,0);

        % outputs
        gm = ['c1',cfg.t1name];
        wm = ['c2',cfg.t1name];
        csf = ['c3',cfg.t1name];

    % ------------------------------------------------------------------------------
    % Generate conservative CSF/WM masks in native space
    % Added following NeuroImage reviews
    % ------------------------------------------------------------------------------
        fprintf(1, '\n\t\t ----- Generating conservative CSF/WM masks ----- \n\n');

        % First, get binary brain mask
        system([cfg.fsldir,'bet ',cfg.t1name,' native_t1_brain -f 0.4 -n -m -R']);
        delete('native_t1_brain.nii')

        % CSF
        % threshold gm and binarise
        system([cfg.fsldir,'fslmaths ',gm,' -thr 0.95 -bin vmask']);
        % dilate twice
        system([cfg.fsldir,'fslmaths vmask -dilD -bin vmask']);
        system([cfg.fsldir,'fslmaths vmask -dilD -bin vmask']);
        % combined with wm and invert
        system([cfg.fsldir,'fslmaths vmask -add ',wm,' -binv vmask']);
        % erode whole brain mask twice
        system([cfg.fsldir,'fslmaths native_t1_brain_mask -eroF e_native_t1_brain_mask']);
        system([cfg.fsldir,'fslmaths e_native_t1_brain_mask -eroF e_native_t1_brain_mask']);
        % multipy eroded brain mask with vmask
        csf = ['v',csf];
        system([cfg.fsldir,'fslmaths vmask -mul e_native_t1_brain_mask ',csf]);
        delete('vmask.nii')
        
        % erode 2 times
        csfname = csf; clear csf;
        csf{1} = [csfname(1:end-4),'_e1.nii'];
        system([cfg.fsldir,'fslmaths ',csfname,' -eroF -bin ',csf{1}]);
        csf{2} = [csfname(1:end-4),'_e2.nii'];
        system([cfg.fsldir,'fslmaths ',csf{1},' -eroF -bin ',csf{2}]);

        % WM
        % mask out non-brain
        system([cfg.fsldir,'fslmaths ',wm,' -mas native_t1_brain_mask b',wm]);
        wm = ['b',wm];
        % erode 5 times
        % we save them all so we can examine how WM to GS correlation drops with each erosion
        % but we only pass the final (5th) erosion to wmmask output variable, which is then fed to noise correction
        wmname = wm; clear wm; wm{1} = [wmname(1:end-4),'_e1.nii'];
        system([cfg.fsldir,'fslmaths ',wmname,' -eroF -bin ',wm{1}]);
        for i = 2:5
            wm{i} = [wmname(1:end-4),'_e',num2str(i),'.nii'];
            system([cfg.fsldir,'fslmaths ',wm{i-1},' -eroF -bin ',wm{i}]);
        end

    % ------------------------------------------------------------------------------
    % Preprocess cfg.EPI
    % ------------------------------------------------------------------------------
        fprintf(1, '\n\t\t ----- Preprocessing EPI image ----- \n\n');
        cd(cfg.preprodir)

    % ------------------------------------------------------------------------------
    % Discard first 4 volumes of cfg.EPI image
    % ------------------------------------------------------------------------------
        % optionally cfg.discard first 4 volumes
        if cfg.discard == 1
            fprintf(1, '\n\t\t ----- Discarding first 4 volumes ----- \n\n');

            outname = ['t',cfg.EPI];
            system([cfg.fsldir,'fslroi ',cfg.rawdir,cfg.EPI,' ',cfg.preprodir,outname,' 4 -1']);
            cfg.EPI = outname;
            tN = cfg.N - 4;
        elseif cfg.discard == 0
            copyfile([cfg.rawdir,cfg.EPI],cfg.preprodir)
            tN = cfg.N;
        end

    % ------------------------------------------------------------------------------
    % Realignment #1
    % Important note: we only retain and use the realignment parameters from this
    % step for use with motion-correction (see Power et al., 2017 PLoS ONE)
    % The realigned data is NOT retained or used further.
    % For that, see below for Realignment #2
    % ------------------------------------------------------------------------------
        fprintf('\n\t\t ----- Realignment ----- \n\n');

        % Make raw_mov dir
        mkdir('raw_mov')
        % Move input file
        movefile(cfg.EPI,'raw_mov')
        cd('raw_mov')

        % Run realignment
        RealignEPI(cfg.EPI,tN)

        % Move input file back
        movefile(cfg.EPI,cfg.preprodir)

        % Delete all in raw_mov dir except rp*.txt
        delete('*.nii','*.mat')

        cd(cfg.preprodir)

    % ------------------------------------------------------------------------------
    % Slice timing correction
    % ------------------------------------------------------------------------------
        % Run optional slice-timing correction
        if cfg.slicetime == 1
            fprintf('\n\t\t ----- Slice-timing correction ----- \n\n');

            % slice-timing correction
            SlicetimeEPI([cfg.preprodir,cfg.EPI], cfg.numSlices, cfg.TR, cfg.order, cfg.refSlice, tN);
        end

    % ------------------------------------------------------------------------------
    % AFNI despike
    % ------------------------------------------------------------------------------
        if cfg.despike == 1
            fprintf('\n\t\t ----- AFNI Despike ----- \n\n');

            % Define inputs to realignment and normalisation
            if cfg.slicetime == 1
                DespikeIn = ['a',cfg.EPI];
                DespikeOut = ['da',cfg.EPI];
            else
                DespikeIn = cfg.EPI;
                DespikeOut = ['d',cfg.EPI];
            end
            
            system([cfg.afnidir,'3dDespike ',DespikeIn]);

            % convert to nifti
            system([cfg.afnidir,'3dAFNItoNIFTI despike+tlrc']);

            % delete afni outputs
            delete('despike+tlrc*')

            % rename output file
            movefile('despike.nii',DespikeOut)
        end

    % ------------------------------------------------------------------------------
    % Realignment #2
    % ------------------------------------------------------------------------------
        fprintf('\n\t\t ----- Realignment ----- \n\n');

        % Define inputs to realignment and normalisation
        if cfg.despike == 1
            RealignIn = DespikeOut;
        elseif cfg.slicetime == 1 & cfg.despike == 0
            RealignIn = ['a',cfg.EPI];
        elseif cfg.slicetime == 0 & cfg.despike == 0
            RealignIn = cfg.EPI;
        end

        RealignEPI(RealignIn,tN)

        RealignOut = ['r',RealignIn];
        meanEPI = ['mean',RealignIn];

    % ------------------------------------------------------------------------------
    % Masking realigned EPI
    % ------------------------------------------------------------------------------
    % In the NYU dataset, I was noticing that the realignment would leave NaNs around
    % edges instead of zeros whenever realignment caused the edges of the image to have
    % no sampled data. These NaNs cause all kinds of issues, so we replace them with zeros
        [hdr,data] = read(RealignOut);
        if any(isnan(data(:))) == 1
            fprintf(1, '\t\t Correcting NaNs in realigned EPI. \n');
            data(isnan(data)) = 0;
            write(hdr,data,RealignOut)
        end

        [hdr,data] = read(meanEPI);
        if any(isnan(data(:))) == 1
            fprintf(1, '\t\t Correcting NaNs in mean EPI. \n');
            data(isnan(data)) = 0;
            write(hdr,data,meanEPI)
        end
        clear hdr data

    % ------------------------------------------------------------------------------
    % Spatial normalisation
        % For some reason, ANTs does not seem to work well with 4D files.
        % So we split the cfg.EPI into 3D files, normalise, concatenate, then clean up.
        % This causes ALOT of junk in the command line...
    % ------------------------------------------------------------------------------
        fprintf('\n\t\t ----- Co-registration and normalisation ----- \n\n');

        % Split 4D file in 3D
        spm_file_split(RealignOut)

        wm_temp = cell(size(wm));
        for i = 1:length(wm_temp)
            wm_temp{i} = [cfg.t1dir,wm{i}];
        end

        csf_temp = cell(size(csf));
        for i = 1:length(csf_temp)
            csf_temp{i} = [cfg.t1dir,csf{i}];
        end

        if isfield(cfg, 't14norm')
            SpatialNormalisationANTs([cfg.preprodir,RealignOut(1:end-4)],tN,[cfg.preprodir,meanEPI],...
                [cfg.t1dir,cfg.t1name],...
                [cfg.t1dir,gm],...
                wm_temp,...
                csf_temp,...
                cfg.mni_template,cfg.antsdir,cfg.funcdir,...
                [cfg.t1dir,cfg.t14norm])
        else ~isfield(cfg, 't14norm')
            SpatialNormalisationANTs([cfg.preprodir,RealignOut(1:end-4)],tN,[cfg.preprodir,meanEPI],...
                [cfg.t1dir,cfg.t1name],...
                [cfg.t1dir,gm],...
                wm_temp,...
                csf_temp,...
                cfg.mni_template,cfg.antsdir,cfg.funcdir)
        end

        % Merge warped cfg.EPI
        system([cfg.fsldir,'fslmerge -t w',RealignOut,' w',RealignOut(1:end-4),'_*.nii']);

        % Clean up 3D files
        delete([RealignOut(1:end-4),'_*.nii'])
        delete(['w',RealignOut(1:end-4),'_*.nii'])

        % EPI outputs
        normEPI = ['w',RealignOut];
        normmeanEPI = ['w',meanEPI];
        
        % T1 outputs
        normT1 = ['w',cfg.t1name];    
        gm = ['w',gm];
        for i = 1:length(wm); wm{i} = ['w',wm{i}]; end
        for i = 1:length(csf); csf{i} = ['w',csf{i}]; end

        % move t1 ouputs to cfg.t1dir
        movefile(normT1,cfg.t1dir)
        movefile(gm,cfg.t1dir)
        for i = 1:length(wm); movefile(wm{i},cfg.t1dir); end
        for i = 1:length(csf); movefile(csf{i},cfg.t1dir); end

    % ------------------------------------------------------------------------------
    % Create binary brain mask
    % ------------------------------------------------------------------------------
        fprintf(1,'\n\t\t ----- Create brain masks ----- \n\n')

        % EPI
        cd(cfg.preprodir)
        MaskIn = normmeanEPI;
        system([cfg.fsldir,'bet ',MaskIn,' epi_brain -f 0.4 -n -m -R']);
        delete('epi_brain.nii')
        epiBrainMask = 'epi_brain_mask.nii';

        % T1
        cd(cfg.t1dir)
        MaskIn = normT1;
        system([cfg.fsldir,'bet ',MaskIn,' t1_brain -f 0.4 -n -m -R']);
        delete('t1_brain.nii')
        t1BrainMask = 't1_brain_mask.nii';

        % Union
        cd(cfg.preprodir)
        system([cfg.fsldir,'fslmaths ',epiBrainMask,' -add ',cfg.t1dir,t1BrainMask,' -bin brain_mask']);
        BrainMask = 'brain_mask.nii';

    % ------------------------------------------------------------------------------
    % Mask out non-brain tissue from EPI image, T1, and GM map
    % ------------------------------------------------------------------------------
        fprintf(1,'\n\t\t ----- Apply brain masks ----- \n\n')
        
        cd(cfg.preprodir)
        % EPI
        system([cfg.fsldir,'fslmaths ',normEPI,' -mas ',cfg.preprodir,BrainMask,' b',normEPI]);
        system([cfg.fsldir,'fslmaths ',normmeanEPI,' -mas ',cfg.preprodir,BrainMask,' b',normmeanEPI]);
        
        % T1
        cd(cfg.t1dir)
        system([cfg.fsldir,'fslmaths ', normT1,' -mas ',cfg.preprodir,BrainMask,' b',normT1]);

        % Tissue masks
        system([cfg.fsldir,'fslmaths ',gm,' -mas ',cfg.preprodir,BrainMask,' b',gm]);
        % system([cfg.fsldir,'fslmaths ',wm,' -mas ',cfg.preprodir,BrainMask,' b',wm]);
        % system([cfg.fsldir,'fslmaths ',csf,' -mas ',cfg.preprodir,BrainMask,' b',csf]);

        % EPI outputs
        normEPI = ['b',normEPI];
        normmeanEPI = ['b',normmeanEPI];
        
        % T1 outputs
        normT1 = ['b',normT1];    
        gm = ['b',gm];

    % ------------------------------------------------------------------------------
    % Generate MNI tissue masks for noise correction and tissue smoothing
    % ------------------------------------------------------------------------------
        cd(cfg.t1dir)

        % Binarise WM/CSF masks following normalisation.
        % Normalisation often introduces non-binary non-zero elements in binary maps
        for i = 1:length(wm); system([cfg.fsldir,'fslmaths ',wm{i},' -bin ',wm{i}]); end
        for i = 1:length(csf); system([cfg.fsldir,'fslmaths ',csf{i},' -bin ',csf{i}]); end

        wmmask = wm{end};
        % Check whether the final eroded WM mask has no voxels
        [hdr,data] = read(wmmask);
        if sum(data(:) > 0) < 5;
            fprintf(1, '\t\t WARNING! the eroded WM mask has too few voxels. Checking previous erosion steps.\n');
            breakloop = 0;
            i = length(wm)-1;
            while breakloop == 0
                [hdr,data_temp] = read(wm{i});
                if sum(data_temp(:) > 0) >= 5;
                    wmmask = wm{i};
                    fprintf(1, '\t\t New WM is erosion %u \n', i);
                    breakloop = 1;
                end
                i = i - 1;
            end
        end

        csfmask = csf{end};
        % Check whether the final eroded CSF mask has no voxels
        [hdr,data] = read(csfmask);
        if sum(data(:) > 0) < 5;
            fprintf(1, '\t\t WARNING! the eroded CSF mask has too few voxels. Checking previous erosion steps.\n');
            [hdr,data_temp] = read(csf{1});
            if sum(data_temp(:) > 0) < 5;
                error('\t\t There are no CSF voxels... please check subject data.\n')
            elseif sum(data_temp(:) > 0) >= 5;
                fprintf(1, '\t\t New CSF is erosion 1 \n');
                csfmask = csf{1};
            end
        end

        gmmask = 'gm50_bin.nii';
        system([cfg.fsldir,'fslmaths ',gm,' -thr .50 -bin ',gmmask]);

        clear hdr data data_temp

    % ------------------------------------------------------------------------------
    % 4D intensity normalisation
    % ------------------------------------------------------------------------------
        IntNormIn = normEPI;

        if cfg.intnorm == 1
            fprintf(1,'\n\t\t ----- EPI intensity normalisation ----- \n\n')
            cd(cfg.preprodir)
            
            IntNormOut = ['i',IntNormIn];
            
            IntensityNormalise(IntNormIn)

            % Get dvars
            dvarsExtract = IntNormOut;
            dvars = GetDVARS(dvarsExtract,BrainMask);
            dlmwrite('dvars.txt',dvars)

        elseif cfg.intnorm == 0
            IntNormOut = IntNormIn;
            dvarsExtract = NaN;
            dvars = NaN;
        end

    % ------------------------------------------------------------------------------
    % Demean & Detrend EPI
    % ------------------------------------------------------------------------------
        if cfg.detr == 1
            fprintf(1,'\n\t\t ----- Detrending ----- \n\n')
            cd(cfg.preprodir)

            DetrendIn = IntNormOut;
            DetrendOut = ['d',DetrendIn];

            % load
            [hdr,data] = read(DetrendIn);
            
            % reshape to 2d
            dim = size(data);
            data = reshape(data,[],dim(4));

            % get the data mean
            if cfg.meanback == 1;
                theMean = mean(data,2);
                theMean = repmat(theMean,[1,tN]);
            end

            % Detrend with all data
            [data_out,~] = JP14_demean_detrend(data);

            % add the mean back
            if cfg.meanback == 1;
                data_out = data_out + theMean;
            end

            data_out = reshape(data_out,dim);
            write(hdr,data_out,DetrendOut)

            % ------------------------------------------------------------------------------
            % Scrubbing (Power et al., 2014. NeuroImage) 
            % ------------------------------------------------------------------------------
            if cfg.intnorm == 1 & cfg.runBandpass == 1
                % get FD Power
                mfile = dir([cfg.preprodir,'raw_mov/rp*.txt']);
                mov = dlmread([cfg.preprodir,'raw_mov/',mfile(1).name]);
                % mfile = dir([cfg.preprodir,'rp*.txt']);
                % mov = dlmread([cfg.preprodir,mfile(1).name]);
                fd = GetFDPower(mov);
                dlmwrite('fdPower.txt',fd)
                
                % Create Power temporal mask
                fdThr = 0.2;
                dvarsThr = 20;
                [scrubmask, exclude] = JP14_GetScrubMask(fd,dvars,cfg.TR,fdThr,dvarsThr);
                dlmwrite('JP14_ScrubMask.txt',scrubmask)

                % Detrend including censor mask
                if exclude == 0
                    [data_out,~] = JP14_demean_detrend(data,~scrubmask(:,2));
                    data_out = reshape(data_out,dim);
                    write(hdr,data_out,['jp14',DetrendOut])
                elseif exclude == 1
                    fprintf(1, '\t\t WARNING! There was not enough uncensored time points to perform scrubbing.\n');
                end
            else
                fdThr = NaN;
                dvarsThr = NaN;
                exclude = NaN;
            end

        elseif cfg.detr == 0
            DetrendOut = IntNormOut;
        end

        clear hdr data data_out

    % ------------------------------------------------------------------------------
    % Spatially smooth the data
    % ------------------------------------------------------------------------------
        fprintf('\n\t\t ----- Spatial smoothing ----- \n\n');
        cd(cfg.preprodir)

        SmoothIn = DetrendOut;
        SmoothEPI(SmoothIn,cfg.kernel,tN)

        if cfg.intnorm == 1 & cfg.runBandpass == 1
            if exclude == 0
                SmoothEPI(['jp14',SmoothIn],cfg.kernel,tN)
            end
        end

    % ------------------------------------------------------------------------------
    % Outputs
    % ------------------------------------------------------------------------------
        outEPI{1} = SmoothIn;
        outEPI{2} = ['s',SmoothIn];

        if cfg.intnorm == 1 & cfg.runBandpass == 1
            if exclude == 0
                outEPI{3} = ['jp14',SmoothIn];
                outEPI{4} = ['sjp14',SmoothIn];
            end
        end

    fprintf('\n\t\t ----- Base preprocessing complete ----- \n\n');
end
