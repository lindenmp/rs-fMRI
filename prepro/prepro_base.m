function [tN,gm,wm,csf,epiBrainMask,t1BrainMask,BrainMask,gmmask,wmmask,csfmask,dvars,dvarsExtract,outEPI] = prepro_base(cfg)
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

    if cfg.detr == 1
        fprintf(1, '\t\t Detrend EPI: yes \n');
    elseif cfg.detr == 0
        fprintf(1, '\t\t Detrend EPI: no \n');
    end

    if cfg.intnorm == 1
        fprintf(1, '\t\t Normalise EPI intensity: yes \n');
    elseif cfg.intnorm == 0
        fprintf(1, '\t\t Normalise EPI intensity: no \n');
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
        fprintf(1, '\n\t\t ----- Segmenting T1 ----- \n\n');
        cd([cfg.datadir,cfg.subject])

        % Clean and reinitialise T1 dir
        % movefile([cfg.t1dir,cfg.t1name],[cfg.datadir,cfg.subject])
        % if exist([cfg.t1dir,'*json']) == 2
        %     movefile([cfg.t1dir,'*json'],[cfg.datadir,cfg.subject])
        % end

        % delete([cfg.t1dir,'*'])

        % movefile([cfg.datadir,cfg.subject,'/',cfg.t1name],cfg.t1dir)
        % if exist([cfg.datadir,cfg.subject,'/*json']) == 2
        %     movefile([cfg.datadir,cfg.subject,'/*json'],cfg.t1dir)
        % end
        
        % cd(cfg.t1dir);

        % First crop out neck
        outname = ['c',cfg.t1name];
        % system([cfg.fsldir,'robustfov -i ',cfg.t1name,' -r ',outname]);
        cfg.t1name = outname;

        % Tissue segment T1 with SPM
        % SegmentT1([cfg.t1dir,cfg.t1name],cfg.spmdir,0,0);

        % outputs
        gm = ['c1',cfg.t1name];
        wm = ['c2',cfg.t1name];
        csf = ['c3',cfg.t1name];

    % ------------------------------------------------------------------------------
    % Preprocess cfg.EPI
    % ------------------------------------------------------------------------------
        fprintf(1, '\n\t\tPreprocessing EPI image\n\n');
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
    % step for use with QCFC benchmarks
    % The realigned data is NOT retained or used further.
    % For that, see below for Realignment #2
    % ------------------------------------------------------------------------------
        fprintf('\n\t\t ----- Realignment ----- \n\n');

        % Make mot dir
        mkdir('mot')
        % Move input file
        movefile(cfg.EPI,'mot')
        cd('mot')

        % Run realignment
        RealignEPI(cfg.EPI,tN)

        % Move input file back
        movefile(cfg.EPI,cfg.preprodir)

        % Delete all in mot dir except rp*.txt
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
    % Spatial normalisation
        % For some reason, ANTs does not seem to work well with 4D files.
        % So we split the cfg.EPI into 3D files, normalise, concatenate, then clean up.
        % This causes ALOT of junk in the command line...
    % ------------------------------------------------------------------------------
        fprintf('\n\t\t ----- Co-registration and normalisation ----- \n\n');

        % Split 4D file in 3D
        spm_file_split(RealignOut)

        SpatialNormalisationANTs([cfg.preprodir,RealignOut(1:end-4)],tN,[cfg.preprodir,meanEPI],...
            [cfg.t1dir,cfg.t1name],...
            [cfg.t1dir,gm],...
            [cfg.t1dir,wm],...
            [cfg.t1dir,csf],...
            cfg.mni_template,cfg.antsdir,cfg.funcdir)

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
        wm = ['w',wm];
        csf = ['w',csf];

        % move t1 ouputs to cfg.t1dir
        movefile(normT1,cfg.t1dir)
        movefile(gm,cfg.t1dir)
        movefile(wm,cfg.t1dir)
        movefile(csf,cfg.t1dir)

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
    % Mask out non-brain tissue from EPI image, T1, and tissue masks
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
        system([cfg.fsldir,'fslmaths ',wm,' -mas ',cfg.preprodir,BrainMask,' b',wm]);
        system([cfg.fsldir,'fslmaths ',csf,' -mas ',cfg.preprodir,BrainMask,' b',csf]);

        % EPI outputs
        normEPI = ['b',normEPI];
        normmeanEPI = ['b',normmeanEPI];
        
        % T1 outputs
        normT1 = ['b',normT1];    
        gm = ['b',gm];
        wm = ['b',wm];
        csf = ['b',csf];

    % ------------------------------------------------------------------------------
    % Generate MNI tissue masks for noise correction and tissue smoothing
    % ------------------------------------------------------------------------------
        cd(cfg.t1dir)

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

        gmmask = 'gm50_bin.nii';
        wmmask = 'wm_final.nii';
        csfmask = 'csf_final.nii';

    % ------------------------------------------------------------------------------
    % Detrend EPI
    % ------------------------------------------------------------------------------
        if cfg.detr == 1
            fprintf(1,'\n\t\t ----- Detrending ----- \n\n')
            cd(cfg.preprodir)

            DetrendIn = normEPI;
            DetrendOut = ['d',DetrendIn];

            % create separate directory for REST detrend function
            dtdir = [cfg.preprodir,'temp'];
            mkdir(dtdir)
            % move 4d file to directory
            movefile([cfg.preprodir,DetrendIn],dtdir,'f')
            
            switch cfg.WhichNii
                case '4D'
                    % detrend epis using rest
                    cd(cfg.preprodir)
                    rest_detrend(dtdir, '_detrend')
                    
                    % Move 4D file back to cfg.preprodir
                    movefile([dtdir,'/',DetrendIn],cfg.preprodir,'f')
                case '3D'
                    % Split 4D file in 3D
                    cd(dtdir)
                    spm_file_split(DetrendIn)
                    
                    % Move 4D file back to cfg.preprodir
                    movefile(DetrendIn,cfg.preprodir,'f')

                    % detrend epis using rest
                    cd(cfg.preprodir)
                    rest_detrend(dtdir, '_detrend')

                    % delete 3D files
                    cd(dtdir)
                    delete([DetrendIn(1:end-4),'*nii'])
            end

            % move output file back to cfg.preprodir
            movefile([dtdir,'_detrend/detrend_4DVolume.nii'],[cfg.preprodir,DetrendOut],'f')

            cd(cfg.preprodir)
            rmdir('temp','s')
            rmdir('temp_detrend','s')
        elseif cfg.detr == 0
            DetrendOut = normEPI;
        end

    % ------------------------------------------------------------------------------
    % 4D intensity normalisation
    % ------------------------------------------------------------------------------
        IntNormIn = DetrendOut;

        if cfg.intnorm == 1
            fprintf(1,'\n\t\t ----- EPI intensity normalisation ----- \n\n')
            cd(cfg.preprodir)
            
            IntNormOut = ['i',IntNormIn];
            
            IntensityNormalise(IntNormIn)

            % Get dvars and save to .mat
            dvarsExtract = IntNormOut;
            dvars = GetDVARS(dvarsExtract,BrainMask);
        elseif cfg.intnorm == 0
            IntNormOut = IntNormIn;
        end

    % ------------------------------------------------------------------------------
    % Spatially smooth the data
    % ------------------------------------------------------------------------------
        fprintf('\n\t\t ----- Spatial smoothing ----- \n\n');
        cd(cfg.preprodir)

        SmoothIn = IntNormOut;

        % Whole brain
        system([cfg.afnidir,'3dBlurInMask -input ',SmoothIn,' -FWHM ',num2str(cfg.kernel),' -mask ',cfg.preprodir,BrainMask,' -prefix s']);
        % convert to nifti
        system([cfg.afnidir,'3dAFNItoNIFTI s+tlrc']);
        % delete afni outputs
        delete('s+tlrc*')
        % rename output file
        movefile('s.nii',['s_',SmoothIn])

        % ------------------------------------------------------------------------------
        % Tissue specific smoothing
        % ------------------------------------------------------------------------------
        % GM
        system([cfg.afnidir,'3dBlurInMask -input ',SmoothIn,' -FWHM ',num2str(cfg.kernel),' -mask ',cfg.t1dir,gmmask,' -prefix sgm']);
        % convert to nifti
        system([cfg.afnidir,'3dAFNItoNIFTI sgm+tlrc']);
        % delete afni outputs
        delete('sgm+tlrc*')
        % rename output file
        movefile('sgm.nii',['sgm_',SmoothIn])       

        % WM
        system([cfg.afnidir,'3dBlurInMask -input ',SmoothIn,' -FWHM ',num2str(cfg.kernel),' -mask ',cfg.t1dir,wmmask,' -prefix swm']);
        % convert to nifti
        system([cfg.afnidir,'3dAFNItoNIFTI swm+tlrc']);
        % delete afni outputs
        delete('swm+tlrc*')
        % rename output file
        movefile('swm.nii',['swm_',SmoothIn])  

        % CSF
        system([cfg.afnidir,'3dBlurInMask -input ',SmoothIn,' -FWHM ',num2str(cfg.kernel),' -mask ',cfg.t1dir,csfmask,' -prefix scsf']);
        % convert to nifti
        system([cfg.afnidir,'3dAFNItoNIFTI scsf+tlrc']);
        % delete afni outputs
        delete('scsf+tlrc*')
        % rename output file
        movefile('scsf.nii',['scsf_',SmoothIn])

    % ------------------------------------------------------------------------------
    % Outputs
    % ------------------------------------------------------------------------------
        outEPI{1} = SmoothIn;
        outEPI{2} = ['s_',SmoothIn];
        outEPI{3} = ['sgm_',SmoothIn];
        outEPI{4} = ['swm_',SmoothIn];
        outEPI{5} = ['scsf_',SmoothIn];

    fprintf('\n\t\t ----- Base preprocessing complete ----- \n\n');
end
