%% IntensityNormalise: 
function [] = IntensityNormalise(rsData)
	% This function performs modal intensity normalisation to the value of 1000
	%
	% ------
	% INPUTS
	% ------
	% rsData 	- location and name of processed resting state ts. e.g., /path/to/dir/epi_brain.nii
	% 			IMPORTANT NOTE, the input .nii MUST have been masked such that non-brain voxels are equal to 0
	% 			e.g., 'fslmaths epi.nii -mas brain_mask.nii epi_brain.nii'
	% 			where brain_mask.nii is generated via something like: 'bet epi.nii brain -f 0.3 -n -m -R'
	% 			see 'Generate brain mask' for more
	% 
	% -------
	% OUTPUTS
	% -------
	% Intensity normalised 4D volume which will be named the same as the input file but with an 'i' prefix
	% to denote intensity normalisation
	%
	% Linden Parkes, Brain & Mental Health Laboratory, 2016
	% ------------------------------------------------------------------------------

	% ------------------------------------------------------------------------------
	% Read in resting state data .nii
	% ------------------------------------------------------------------------------
	[hdr,ts] = read(rsData);

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
	% Generate brain mask
	% ------------------------------------------------------------------------------
	% I do this by finding voxels (columns) that are all-zero
	% This is why it is important to input 4d epi file that has brain masked BEFORE
	% inputting to this function
	% Leaving these voxels in will guarantee a mode value of 0...not useful...
	
	% Find zeros voxels
	mask = ts == 0;
	% sum across rows (time points) to find voxels with all-zero voxels
	mask = sum(mask);
	% Retain only those sums that equal length of times
	mask(mask == N) = 1;
	% Invert to brain voxels
	mask = ~mask;
	% Make logical
	mask = logical(mask);

	% ------------------------------------------------------------------------------
	% Mask out non brain voxels
	% ------------------------------------------------------------------------------
	ts = ts(:,mask);

	% ------------------------------------------------------------------------------
	% Get modal value
	% ------------------------------------------------------------------------------
	tsMode = mode(ts(:));

	% ------------------------------------------------------------------------------
	% Normalise by mode
	% ------------------------------------------------------------------------------
	if tsMode == 0
		fprintf(1, 'NOTE, modal value is 0! exiting...\n');
	else
		ts = (ts./tsMode)*1000;
	end

	% ------------------------------------------------------------------------------
	% Reshape back to 4D
	% ------------------------------------------------------------------------------
	% Add non brain voxels back in
	numVoxels = numel(mask);
	tsTemp = zeros(N,numVoxels);
	
	idx = find(mask);
	numBrainVoxels = length(idx);

	% loop over each brain voxel column
	for i = 1:numBrainVoxels
		% write it to the corresponding column in tsTemp
		tsTemp(:,idx(i)) = ts(:,i);
	end

	tsOut = tsTemp;

	% ------------------------------------------------------------------------------
	% Write out filtered ts
	% ------------------------------------------------------------------------------
	% Put time back on second dimensions
	tsOut = tsOut';
	% Reshape back to 4D matrix
	tsOut = reshape(tsOut,dim(1),dim(2),dim(3),dim(4));

	% Write
	[fPath,fName,fExt] = fileparts(rsData);
	outName = ['i',fName,fExt];
	write(hdr,tsOut,outName)

end