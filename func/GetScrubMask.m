%% GetScrubMask: 
function Mask = GetScrubMask(fd,dvars,fdThr,dvarsThr,scrubProximal)

	% This function generates a mask for use with scrubbing/censoring as in Power's work
	% Mask is generated based on concurrent suprathreshold changes in both FD & dvars
	% N-1, N+1, and N+2 volumes are also flagged surrounding these 'bad' volumes
	%
	% ------
	% INPUTS
	% ------
	% fd 		- fd metric from Power's work (see GetFDPower.m)
	% dvars 	- see GetDVARS.m. Note, should run on only minimally preprocessed data (e.g., before noise correction)
	% 
	% fdThr 	- Threshold for detecting FD movements. default = 0.2
	% dvarsThr 	- Threshold for detecting DVARS changes. default = 0.3
	%
	% -------
	% OUTPUTS
	% -------
	% Mask  	- Logical vector for use with censoring/scrubbing of fully preprocessed data
	%
	% Linden Parkes, Brain & Mental Health Laboratory, 2016
	% ------------------------------------------------------------------------------

	if nargin < 2
		fdThr = 0.2;
	end

	if nargin < 3
		dvarsThr = 3;
	end

	if nargin < 4;
		scrubProximal = 'yes';
	end
	
	N1 = length(fd);
	N2 = length(dvars);
	if N1 ~= N2
		fprintf(1, 'WARNING! fd and dvars are different lengths');
	else
		N = N1;
	end
	
	% Genereate FD mask
	fdMask = fd > fdThr;

	% Genereate dvars mask
	dvarsMask = dvars > dvarsThr;

	% Combine
	MaskTemp = fdMask + dvarsMask;

	% Find volumes = 2. i.e., when both FD and dvars are suprathreshold
	Mask = MaskTemp >= 2;

	switch scrubProximal
		case 'yes'
			% Also flag 1 back and 2 forward (as in Power 2012)
			idx = find(Mask);
			for i = idx
				% 1 back
				Mask(i-1) = 1;

				% 2 forward
				Mask(i+1) = 1;
				Mask(i+2) = 1;
			end

			% Check if there is any overlap
			if any(Mask >= 2)
				idx = find(Mask >= 2);
				Mask(idx) = 1;
			end

			% Change to logical
			Mask = logical(Mask);

			% Check if length of mask exceeds numVols
			% This happens where volumes near the end of acquisition are marked for scrubbing
			if length(Mask) > N
				% Trim off extra mask volumes
				Mask(N+1:end) = [];
			end
	end

end