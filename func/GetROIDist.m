%% GetROIDist: 
function [ROIDist] = GetROIDist(ROICoords)
	% This function produces pairwise distances between a list of ROIs in MNI space
	%
	% ------
	% INPUTS
	% ------
	% ROICoords 	- a numROIs x 3 (xyz) matrix of MNI coordinates:
	% 				column 1 = x
	% 				column 2 = y
	% 				column 3 = z
	%
	% -------
	% OUTPUTS
	% -------
	% ROIDist 		- a numROIs x numROIs matrix of pairwise distances between ROIs
	%
	% Linden Parkes, Brain & Mental Health Laboratory, 2016
	% ------------------------------------------------------------------------------
	numROIs = size(ROICoords,1);

    for i = 1:numROIs
    	for j = 1:numROIs
    		x1 = ROICoords(i,1);
    		y1 = ROICoords(i,2);
    		z1 = ROICoords(i,3);

    		x2 = ROICoords(j,1);
    		y2 = ROICoords(j,2);
    		z2 = ROICoords(j,3);

    		ROIDist(i,j) = sqrt((x1-x2)^2 + (y1-y2)^2 + (z1-z2)^2);
    	end
    end

end