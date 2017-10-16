function [ts_compartment,key_compartment] = GetTSCompartment(fsldir,rsData,gmMask,wmMask)
%% GetTSCompartment
	% This function will output a 2d time series matrix (Voxels x Time) for processed resting state data
	% Rows of the matrix are grouped by tissue compartments (grey matter, white matter, csf)
	% Tissue probability masks from SPM8's segmentation toolbox are required

	% NOTE. Dimensions of resting state data and tissue maps must be the same

	% ------
	% INPUTS
	% ------
	% 
	% fsldir - location of fsl binaries. e.g., fsldir = '/usr/share/fsl/5.0/bin/';
	% rsData - location and name of processed resting state data. e.g., /path/to/dir/epi_prepro.nii
	% gmMask - location and name of grey matter probability mask. e.g., /path/to/dir/crwc1t1.nii
	% wmMask - location and name of white matter probability mask. e.g., /path/to/dir/crwc2t1.nii
	% maskThr - (optional) value between 0 and 1 that is used to threshold and binarise the tissue maps. default = 0.50
	% 
	% -------
	% OUTPUTS
	% -------
	% 
	% ts_compartment - 2d time series matrix (Voxels x Time)
	% key_compartment - vector denoting which compartment a given voxel belongs to. may be used for plotting
		% 1 = gm
		% 2 = wm

	% if nargin < 6
	% 	maskThr = 0.50;
	% end

	% ------------------------------------------------------------------------------
	% set FSL environments 
	% ------------------------------------------------------------------------------
	setenv('FSLDIR',fsldir(1:end-4));
	setenv('FSLOUTPUTTYPE','NIFTI');

	% setenv('LD_LIBRARY_PATH',[getenv('PATH'),getenv('LD_LIBRARY_PATH'),':/usr/lib/fsl/5.0'])
	% setenv('/usr/lib/fsl/5.0')

	% ------------------------------------------------------------------------------
	% Read in resting state data .nii
	% ------------------------------------------------------------------------------
	[~,ts] = read(rsData);

	% ------------------------------------------------------------------------------
	% Load in tissue masks
	% ------------------------------------------------------------------------------
	[~,gm] = read(gmMask);
	[~,wm] = read(wmMask);

	% Binarise just incase binary masks weren't input
	gm(gm > 0) = 1;
	wm(wm > 0) = 1;

	% % ------------------------------------------------------------------------------
	% % Threshold tissue masks
	% % ------------------------------------------------------------------------------
	% system([fsldir,'fslmaths ',gmMask,' -thr ',num2str(maskThr),' -bin gm_temp']);
	% system([fsldir,'fslmaths ',wmMask,' -thr ',num2str(maskThr),' -bin wm_temp']);

	% % ------------------------------------------------------------------------------
	% % Load in thresholded tissue masks
	% % ------------------------------------------------------------------------------
	% [~,gm] = read('gm_temp.nii');
	% [~,wm] = read('wm_temp.nii');

	% delete('gm_temp.nii','wm_temp.nii')

	% ------------------------------------------------------------------------------
	% Reshape in 2d time series matrix
	% ------------------------------------------------------------------------------
	dim = size(ts);
    T = dim(4);
	ts = reshape(ts,[],T);
	% Put time on first dimension
	ts = ts';
	
	% ------------------------------------------------------------------------------
	% Reshape tissue masks into 2d vectors
	% ------------------------------------------------------------------------------
	% tissue masks
	gm = reshape(gm,[],1);
	gm = gm'; gm = logical(gm);

	wm = reshape(wm,[],1);
	wm = wm'; wm = logical(wm);

	% ------------------------------------------------------------------------------
	% Remove overlap between tissue masks to avoid duplicate voxels in ts
	% ------------------------------------------------------------------------------
	% create overlap vector
	overlap = [gm; wm];
	overlap = sum(overlap);

	if numel(unique(overlap)) > 2
		overlap = overlap > 1;

		gm(overlap) = 0;
		wm(overlap) = 0;
	end

	% ------------------------------------------------------------------------------
	% Generate 2d ts matrix orgasinised by tissue compartments
	% ------------------------------------------------------------------------------
	ts_compartment = [ts(:,gm) ts(:,wm)];

	% create index variable that denotes tissue compartment of each row in ts_comparment
	key_compartment = [ones(1,sum(gm == 1)) - 0.5, ones(1,sum(wm == 1))];

end