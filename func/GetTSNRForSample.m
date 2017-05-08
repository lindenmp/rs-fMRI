%% GetTSNRForSample: 
function [tSNR,tSNR_ROI] = GetTSNRForSample(cfg,numSubs,numROIs,Parc)
	% ------------------------------------------------------------------------------
	% 
	% Linden Parkes, Brain & Mental Health Laboratory, 2016
	% ------------------------------------------------------------------------------

	% ------------------------------------------------------------------------------
	% Containers
	% ------------------------------------------------------------------------------
	tSNR = zeros(numSubs,1);

	tSNR_ROI = zeros(numSubs,numROIs);

	% ------------------------------------------------------------------------------
	% Loop
	% ------------------------------------------------------------------------------
	for i = 1:numSubs
		% msg = sprintf('\tProcessing subject %u/%u: %s\n',i,numSubs,ParticipantIDs{i});
	    % fprintf(msg);

	    % get ts
	    ts = cfg(i).roiTS{Parc};

	    % ensure timeseries are non-negative
	    constant = 1000;
	    ts = ts + constant;
	    if any(ts(:) < 0)
	    	error('There a negative values in time series data\n')
	    end

	    % time series mean of each ROI
	    tMean = mean(ts,1);
	    % time series standard deviation of each ROI
	    tSd = std(ts,1);
	    % ratio of mean to sd for each ROI
	    tSNR_ROI(i,:) = tMean./tSd;
	    % average over ROIs
	    tSNR(i) = mean(tSNR_ROI(i,:));

	    % n = numel(msg);
	    % fprintf(repmat('\b',1,n));
	end
end