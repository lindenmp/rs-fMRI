%% GetScrubMask: 
function [ScrubMask, exclude] = JP12_GetScrubMask(fd,dvars,TR,fdThr,dvarsThr,scrubProximal)

	% This function generates a mask for use with scrubbing/censoring as in Power's work
	% ScrubMask is generated based on concurrent suprathreshold changes in both FD & dvars
	% N-1, N+1, and N+2 volumes are also flagged surrounding these 'bad' volumes
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
	% -------
	% OUTPUTS
	% -------
	% ScrubMask - Logical vector for use with censoring/scrubbing of fully preprocessed data
	%
	% Linden Parkes, Brain & Mental Health Laboratory, 2016
	% ------------------------------------------------------------------------------

	if nargin < 4
		fdThr = 0.2;
	end

	if nargin < 5
		dvarsThr = 30;
	end

	if nargin < 6;
		scrubProximal = 'no';
	end
	
	N1 = length(fd);
	N2 = length(dvars);
	if N1 ~= N2
		fprintf(1, 'WARNING! fd and dvars are different lengths');
	else
		numVols = N1;
	end
	
	% Genereate FD ScrubMask
	fdMask = fd > fdThr;

	% Genereate dvars ScrubMask
	dvarsMask = dvars > dvarsThr;

	% Combine
	tempMask = fdMask + dvarsMask;

	% Find volumes = 2. i.e., when FD OR dvars are suprathreshold (Power 2013)
	ScrubMask = tempMask >= 1;

	switch scrubProximal
		case 'yes'
			% Also flag 1 back and 2 forward (as in Power 2012)
			idx = find(ScrubMask);
			for i = idx
				% 1 back
				ScrubMask(i-1) = 1;

				% 2 forward
				ScrubMask(i+1) = 1;
				ScrubMask(i+2) = 1;
			end

			% Check if there is any overlap
			if any(ScrubMask >= 2)
				idx = find(ScrubMask >= 2);
				ScrubMask(idx) = 1;
			end

			% Change to logical
			ScrubMask = logical(ScrubMask);

			% Check if length of ScrubMask exceeds numVols
			% This happens where volumes near the end of acquisition are marked for scrubbing
			if length(ScrubMask) > numVols
				% Trim off extra ScrubMask volumes
				ScrubMask(numVols+1:end) = [];
			end
	end

	% ------------------------------------------------------------------------------
	% Determine exclusion
	% ------------------------------------------------------------------------------
	% threshold for exclusion in minutes
	thresh = 4;
	
	% number of volumes - number of scrubbed volumes
	numCVols = numVols - sum(ScrubMask);

	NTime = (numCVols * TR)/60;

	if NTime < thresh
		exclude = 1;
	else
		exclude = 0;
	end

	exclude = logical(exclude);

end
