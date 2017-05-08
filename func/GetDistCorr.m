%% GetDistCorr: 
function [QC_FC,P] = GetDistCorr(QC,FC)
	% This function correlates functional connectivity matrices with some measure of movement/motion across subjects
	%
	%
	% ------
	% INPUTS
	% ------
	% QC			- a 1xN vector containing a summary metric of movement for each of N subjects. 
	% 				e.g., Yan et al. (2013) summary of framewise displacement
	% FC 			- a numROIs x numROIs x numSubs matrix containing functional connecivity matrices for each subject
	% 				NOTE, the third dimension of the FC matrix must match the QC vector
	%
	% -------
	% OUTPUTS
	% -------
	% QC_FC 		- a numROIs x numROIs matrix of correlations between QC and FC across subjects
	% P 			- a numROIs x numROIs matrix of the significance of the above matrix
	%
	% Linden Parkes, Brain & Mental Health Laboratory, 2016
	% ------------------------------------------------------------------------------

	numROIs = size(FC,1);
	numSubs = size(FC,3);

	% reshape into 2d matrix with subjects on first dimension
	FC = reshape(FC,numROIs*numROIs,numSubs);
	FC = FC';

	% correlate QC with FC
	[QC_FC,P] = corr(QC,FC);

	% reshape to 2d matrix
	QC_FC = reshape(QC_FC,numROIs,numROIs);
	P = reshape(P,numROIs,numROIs);

end