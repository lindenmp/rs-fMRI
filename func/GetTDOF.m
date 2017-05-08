%% GetTDOF: 
function [tDOF_mean,tDOF_std] = GetTDOF(cfg,datadir,ParticipantIDs,preprostr,removeNoise)
	% ------------------------------------------------------------------------------
	% Get tDOF
	% ------------------------------------------------------------------------------
	fprintf(1, 'Computing tDOF: %s\n',removeNoise);
	tDOFtemp = zeros(numSubs,1);
	for j = 1:numSubs

		% get tDOF
		% First, find size of second dimension of noiseTS
		if runSR
			tDOFtemp(j) = size(cfg(j).noiseTS_spikereg,2);
		else
			tDOFtemp(j) = size(cfg(j).noiseTS,2);
		end

		if runScrub
			tDOFtemp(j) = tDOFtemp(j) + sum(ScrubMask{j});
		end

		% Then, if ICA-AROMA pipeline, find number of ICs and add to tDOF
		if ~isempty(strfind(removeNoise,'ICA-AROMA'))
			if runSR
				x = dlmread([datadir,ParticipantIDs{j},preprostr,removeNoise,'+SpikeReg/classified_motion_ICs.txt']);
			else
				x = dlmread([datadir,ParticipantIDs{j},preprostr,removeNoise,'/classified_motion_ICs.txt']);
			end
			tDOFtemp(j) = tDOFtemp(j) + length(x);
		end
	end

	% ------------------------------------------------------------------------------
	% Calculate mean temporal degrees of freedom lost
	% ------------------------------------------------------------------------------
	% tDOF will be the same for most pipelines, but some have variable regressor amounts
	% So we take mean over subjects
	tDOF_mean(i) = mean(tDOFtemp);
	tDOF_std(i) = std(tDOFtemp);