function [dvars] = GetDVARS(rsData,MaskIn)

	% This function computes derivate variance across voxels as detailed in 
	% Power et al. (2012) NeuroImage and Power et al. (2014) NeuroImage.
	%
	% ------
	% INPUTS
	% ------
	% rsData - location and name of 4D resting state data. e.g., /path/to/dir/epi.nii
	%
	% MaskIn - location and name of binary brain mask . e.g., /path/to/dir/brain_mask.nii
	% 
	% -------
	% OUTPUTS
	% -------
	% dvars - Derivative VARiance.
	%
	% Linden Parkes, Brain & Mental Health Laboratory, 2016
	% ------------------------------------------------------------------------------
    
	% ------------------------------------------------------------------------------
	% Read in resting state data .nii
	% ------------------------------------------------------------------------------
	[~,ts] = read(rsData);

	% ------------------------------------------------------------------------------
	% Reshape in 2d time series matrix
	% ------------------------------------------------------------------------------
	% rsData
	dim = size(ts);
    % number of time points
    N = dim(4);
	ts = reshape(ts,[],N);
	% Put time on first dimension
	ts = ts';

	% ------------------------------------------------------------------------------
	% Load in brain mask
	% ------------------------------------------------------------------------------
	[hdr_mask,mask] = read(MaskIn);

	% reshape
	mask = reshape(mask,[],1);
	mask = mask';

	% binarise if not already
	if numel(unique(mask)) > 2
		mask(mask > 0) = 1;
	end
	
	mask = logical(mask);

	% ------------------------------------------------------------------------------
	% Mask out non brain voxels
	% ------------------------------------------------------------------------------
	ts = ts(:,mask);

	% ------------------------------------------------------------------------------
	% Compute backwards differences
	% ------------------------------------------------------------------------------
	ts = [
		zeros(1,size(ts,2));
		diff(ts)
		];

	% Root Mean Square across voxels
    dvars = rms(ts,2);

end