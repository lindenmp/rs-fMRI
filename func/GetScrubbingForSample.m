%% GetScrubbingForSample: 
function [mask,mov,fdPower,dvars] = GetScrubbingForSample(datadir,ParticipantIDs,str)
	% ------------------------------------------------------------------------------
	% This script uses the movement parameters generate during SPM8's realignment
	% to calculate framewise displacement according to perform scrubbing as in Power's work
	% 
	% Linden Parkes, Brain & Mental Health Laboratory, 2016
	% ------------------------------------------------------------------------------

	numSubs = length(ParticipantIDs);
	% ------------------------------------------------------------------------------
	% Containers
	% ------------------------------------------------------------------------------
	mask = cell(numSubs,1);

	mov = cell(numSubs,1);

	fdPower = cell(numSubs,1);

	dvars = cell(numSubs,1);

	% ------------------------------------------------------------------------------
	% Loop
	% ------------------------------------------------------------------------------
	for i = 1:numSubs
		% msg = sprintf('\tProcessing subject %u/%u: %s\n',i,numSubs,ParticipantIDs{i});
	    % fprintf(msg);

		workdir = [datadir,ParticipantIDs{i},str];

		% ------------------------------------------------------------------------------
		% Compute fd (Power2012)
		% ------------------------------------------------------------------------------
	    % read in motion
	    mfile = dir([workdir,'rp*txt']);
	    mov{i} = dlmread([workdir,mfile(1).name]);
		
		% Threshold for flagging problem volumes
		fdPowerThr = 0.2;

		% Get FD
		fdPower{i} = GetFDPower(mov{i});

		% ------------------------------------------------------------------------------
		% Compute DVARS
		% NOTE: Important that the epi used for this is the one BEFORE major noise correction
		% Here I use the intensity normalised detrended epi output from core image preprocessing
		% ------------------------------------------------------------------------------
	    dvarsThr = 20;

	    % Get dvars
	    dvars{i} = dlmread([workdir,'dvars.txt']);

		% ------------------------------------------------------------------------------
		% Generate scrubbing mask
		% ------------------------------------------------------------------------------
		% scrubProximal = 'yes';
		scrubProximal = 'no';
		mask{i} = JP12_GetScrubMask(fdPower{i},dvars{i},fdPowerThr,dvarsThr,scrubProximal);
	    
	    % n = numel(msg);
	    % fprintf(repmat('\b',1,n));
	end

end