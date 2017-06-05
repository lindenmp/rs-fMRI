%% GetSpikeRegressors: 
function [spikereg] = GetSpikeRegressors(fd,thr)

	if nargin < 2
		thr = 0.25;
	end

	% find indices of spikes
	if any(fd > thr)
	    % fprintf(1,'Generating spike regressors...');

		idx = find(fd > thr);

		% find length of time series.
		N = length(fd);

		% find number of spikes
		Nspike = numel(idx);

		% fprintf(1, 'found %u spikes...',Nspike);

		% initialise spike regression matrix
		spikereg = zeros(N,Nspike);

		for i = 1:Nspike
			spikereg(idx(i),i) = 1;
		end
	
		% fprintf(1,' done\n');
	else
	    % fprintf(1,'no spikes detected\n');
	    spikereg = [];
	end
end
