%% GetScrubMask: 
function [ScrubMask, exclude] = GetScrubMask(fd,dvars,TR,fdThr,dvarsThr)

	% This function generates a ScrubMask for use with scrubbing/censoring as in Power's 2014 work
	% ScrubMask is generated based on concurrent suprathreshold changes in both FD & dvars
	%
	% ------
	% INPUTS
	% ------
	% fd 		- fd metric from Power's work (see GetFDPower.m)
	% dvars 	- see GetDVARS.m. Note, should run on only minimally preprocessed data (e.g., before noise correction)
	% TR 		- TR in seconds
	% 
	% fdThr 	- Threshold for detecting FD movements. default = 0.2
	% dvarsThr 	- Threshold for detecting DVARS changes. default = 20
	% 
	%
	% -------
	% OUTPUTS
	% -------
	% ScrubMask  	- Logical vector for use with censoring/scrubbing of fully preprocessed data
	% exclude 	- a binary variable where 1 = exclude subject, 0 = retain
	%
	% Linden Parkes, Brain & Mental Health Laboratory, 2017
	% ------------------------------------------------------------------------------

	if nargin < 4
		fdThr = 0.2;
	end

	if nargin < 5
		dvarsThr = 20;
	end

	N1 = length(fd);
	N2 = length(dvars);
	if N1 ~= N2
		fprintf(1, 'WARNING! fd and dvars are different lengths');
	else
		numVols = N1;
	end
	
	% ------------------------------------------------------------------------------
	% Step 1: identify suprathreshold FD and DVARS
	% ------------------------------------------------------------------------------
	% Genereate FD ScrubMask
	fdMask = fd > fdThr;

	% Genereate dvars ScrubMask
	dvarsMask = dvars > dvarsThr;

	% Combine
	tempMask = fdMask + dvarsMask;

	% Find volumes = 1. i.e., when FD OR dvars are suprathreshold
	ScrubMask = tempMask >= 1;

	% ------------------------------------------------------------------------------
	% Step 2: identify uncensored segmenets with < 5 contiguous volumes
	% ------------------------------------------------------------------------------
	paddedMask = [1; ScrubMask(:,1); 1];
	idx = find(paddedMask);
	idx = idx - 1; idx(1) = 1;

	tempMask = zeros(size(fd));
	for i = 1:length(idx)-1
		start_idx = idx(i);
		end_idx = idx(i+1);
		if length(start_idx:end_idx)-2 < 5
			tempMask(start_idx:end_idx-1,1) = 1;
		end
	end

	% add overlap 
	tempMask = tempMask + ScrubMask(:,1);
	tempMask(tempMask ~= 0) = 1;

	ScrubMask(:,2) = tempMask;

	% ------------------------------------------------------------------------------
	% Step 3: identify runs with < 50 remaining volumes
	% Power's work included datasets with multiple scans (runs). We only have a
	% single run per subject (except NYU)
	% Not relevant to us since we do not have more than one run in our primary datasets
	% ------------------------------------------------------------------------------

	% ------------------------------------------------------------------------------
	% Step 4: identify subjects with <4 mins uncensored volumes across runs
	% Note, we only have a single run
	% Note, we use a more lenient threshold of 4-minutes instead of 5-minutes as per Satterthwaite et al. 2013
	% ------------------------------------------------------------------------------
	% threshold for exclusion in minutes
	thresh = 4;

	numCVols = numVols - sum(ScrubMask(:,2));

	NTime = (numCVols * TR)/60;

	if NTime < thresh
		exclude = 1;
	else
		exclude = 0;
	end

	exclude = logical(exclude);
end
