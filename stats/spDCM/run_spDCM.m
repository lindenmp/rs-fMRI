function [] = run_spDCM(datadir,data,units,N,TR,TE,roidir)
	% This script is designed to replicate the DCM for resting state tutorial in the
	% SPM12 manual.
	% I left out the parts where noise time series are extracted and then used to 
	% adjust the SPM. Thus, the data input to this function has to have been cleaned
	% already.

	% Copyright (C) 2017, Linden Parkes <lindenparkes@gmail.com>,

	% ------------------------------------------------------------------------------
	% Setup output directory
	% ------------------------------------------------------------------------------
	dcmdir = [datadir,'spDCM/'];
	if exist(dcmdir) == 0
		fprintf(1,'Initialising dcmdir\n')
		mkdir(dcmdir)
	elseif exist(dcmdir) == 7
		fprintf(1,'Cleaning and re-initialising dcmdir\n')
		rmdir(dcmdir,'s')
		mkdir(dcmdir)
	end

	% ------------------------------------------------------------------------------
	% Initiate SPM
	% ------------------------------------------------------------------------------
	spm('defaults','fmri');
	spm_jobman('initcfg')
	spm_get_defaults('mask.thresh',-Inf)

	%-----------------------------------------------------------------------
	% Estimate GLM
	%-----------------------------------------------------------------------

	matlabbatch{1}.spm.stats.fmri_spec.dir = {dcmdir};
	matlabbatch{1}.spm.stats.fmri_spec.timing.units = units;
	matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TR;
	matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
	matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
	%%

	for i = 1:N
	    matlabbatch{1}.spm.stats.fmri_spec.sess.scans{i,1} = [datadir,data,',',num2str(i)];
	end

	%%
	matlabbatch{1}.spm.stats.fmri_spec.sess.cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
	matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {''};
	matlabbatch{1}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
	matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {''};
	matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = 128;
	matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
	matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
	matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
	matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
	matlabbatch{1}.spm.stats.fmri_spec.mthresh = -Inf;
	matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
	matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
	matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
	matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
	matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

	spm_jobman('run',matlabbatch);
	clear matlabbatch

	%-----------------------------------------------------------------------
	% Generate VOIs
	%-----------------------------------------------------------------------
	ROIs = dir([roidir,'*.nii']);
	numROIs = size(ROIs,1);

	for i = 1:numROIs

		matlabbatch{1}.spm.util.voi.spmmat = {[dcmdir,'SPM.mat']};
		matlabbatch{1}.spm.util.voi.adjust = NaN;
		matlabbatch{1}.spm.util.voi.session = 1;
		matlabbatch{1}.spm.util.voi.name = ROIs(i).name(1:end-4);
		matlabbatch{1}.spm.util.voi.roi{1}.mask.image = {[roidir,ROIs(i).name,',1']};
		matlabbatch{1}.spm.util.voi.roi{1}.mask.threshold = 0;
		matlabbatch{1}.spm.util.voi.roi{2}.mask.image = {[dcmdir,'mask.nii,1']};
		matlabbatch{1}.spm.util.voi.roi{2}.mask.threshold = 0;
		matlabbatch{1}.spm.util.voi.expression = 'i1&i2';

		spm_jobman('run',matlabbatch);
		clear matlabbatch
	end

	% ------------------------------------------------------------------------------
	% Specify DCM
	% There wasn't any obvious way to do this with the SPM command line tools that 
	% didn't require interaction with the SPM GUI (e.g., spm_dcm_specify) and the SPM12
	% manual only gives instructions for specifying via the GUI not via the batch editor.
	% 
	% So, I made a template DCM.mat file using the GUI and then wrote the code below
	% to replicate that template. This way the code can run without being dependent 
	% on a .mat file
	% ------------------------------------------------------------------------------
	load([dcmdir,'SPM.mat'])

	clear DCM

	% Load regions of interest
	for i = 1:numROIs
		load(fullfile(dcmdir,['VOI_',ROIs(i).name(1:end-4),'_1.mat']),'xY');
		DCM.xY(i) = xY;
	end
	clear xY

	DCM.n = length(DCM.xY); % number of regions. I could have pulled this from numROIs
	DCM.v = length(DCM.xY(1).u); % number of time points. Could also pull this from function input N

	% Define connections
	DCM.a = ones(DCM.n,DCM.n); % this builds fully connected model
	DCM.b = [];
	DCM.c = [];
	DCM.d = [];

	% Experimental inputs (which there are none)
	DCM.U.u = zeros(DCM.v,1);
	DCM.U.name = {'null'};

	% Time series
	DCM.Y.dt  = SPM.xY.RT;
	DCM.Y.X0  = DCM.xY(1).X0;
	for i = 1:DCM.n
	    DCM.Y.y(:,i)  = DCM.xY(i).u;
	    DCM.Y.name{i} = DCM.xY(i).name;
	end

	DCM.Y.Q    = spm_Ce(ones(1,DCM.n)*DCM.v);

	% DCM parameters and options
	DCM.delays = repmat(SPM.xY.RT/2,DCM.n,1);
	DCM.TE     = TE;

	DCM.options.nonlinear = 0;
	DCM.options.two_state = 0;
	DCM.options.stochastic = 1;
	DCM.options.centre = 1;

	% Save
	save([dcmdir,'DCM.mat'],'DCM')

	% ------------------------------------------------------------------------------
	% Estimate DCM
	% ------------------------------------------------------------------------------
	matlabbatch{1}.spm.dcm.fmri.estimate.dcmmat = {[dcmdir,'DCM.mat']};
	matlabbatch{1}.spm.dcm.fmri.estimate.analysis = 'csd';

	spm_jobman('run',matlabbatch);
	clear matlabbatch

end

