%% GetExcludeForSample: 
function [exclude,mov,fdJenk,fdJenk_m,fdPower,fdPower_m,dvars,JP12ScrubMask,JP14ScrubMask] = GetExcludeForSample(datadir,ParticipantIDs,TR,str,mname)
	% ------------------------------------------------------------------------------
	% This script uses the movement parameters generate during SPM8's realignment
	% to calculate framewise displacement according to two methods, (1) Power and
	% (2) Jenkinson. Then we apply the following criteria to determine to Jenkinson's
	% FD measure to determine whether a participant ought to be excluded:
	% 1) if mean (over time) FD exceeds 0.2mm
	% OR
	% 2) if greater than 20% of FDs exceed a threshold of 0.25mm
	% OR
	% 3) if any single FD is greater than 5mm
	% If any of the above criteria are satisfied, a subject is marked for exclusion
	% 
	% 
	% INPUT
	% datadir 	- string containing directory to where subject's folder are located
    %   			datadir = '/gpfs/M2Home/projects/Monash076/Linden/OCDPG/data/';
    % ParticipantIDs 	- cell containing subject identifiers
	% str 		- string containing subfolder structure that points to where the realignment params are
	%				str = '/rfMRI/prepro/';
	% 				if the realign params are immediately inside the subject's directory, then make str = '';
	% 
	% All three of these input strings should concatenate to form the path straight for a subject's realignment params
	% see 'workdir' below
	% workdir = [datadir,ParticipantIDs{i},str];
	% 
	% PRIMARY OUTPUT
	% exclude - a logical vector denoting whether a subject should be excluded
	% 
	% Linden Parkes, Brain & Mental Health Laboratory, 2016
	% ------------------------------------------------------------------------------

	if nargin < 4
		str = '';
	end

	if nargin < 5
		mname = 'rp*txt';
	end

	numSubs = length(ParticipantIDs);

	% ------------------------------------------------------------------------------
	% Containers
	% ------------------------------------------------------------------------------
	mov = cell(numSubs,1);

	fdJenk = cell(numSubs,1);
	fdJenk_m = zeros(numSubs,1);

	fdPower = cell(numSubs,1);
	fdPower_m = zeros(numSubs,1);

	dvars = cell(numSubs,1);

	JP12ScrubMask = cell(numSubs,1);
	JP14ScrubMask = cell(numSubs,1);

	exclude = zeros(numSubs,5);

	% ------------------------------------------------------------------------------
	% Loop
	% ------------------------------------------------------------------------------
	for i = 1:numSubs
		msg = sprintf('\tProcessing subject %u/%u: %s\n',i,numSubs,ParticipantIDs{i});
	    fprintf(msg);

	    % ------------------------------------------------------------------------------
		% Load in movement parameters from realignment
		% ------------------------------------------------------------------------------
		workdir = [datadir,ParticipantIDs{i},str];
	    mfile = dir([workdir,'raw_mov/',mname]);
	    mov{i} = dlmread([workdir,'raw_mov/',mfile(1).name]);

	    [pathstr,name,ext] = fileparts(mfile.name);
	    if ismember('.par',ext,'rows')
	    	% If movement file has .par extension, assume FSL organisation of columns (rot,trans) and reordered to SPM organisation (trans,rot)
	    	mov{i} = mov{i}(:,[4:6,1:3]);
	    end

		numVols = size(mov{i},1);

	    % ------------------------------------------------------------------------------
		% Compute fd Jenkinson (used by Satterthwaite2012)
	    % ------------------------------------------------------------------------------
		% Threshold for detecting 'spikes'
		fdJenkThr = 0.25;

		% Get FD
		fdJenk{i} = GetFDJenk(mov{i});

		% Calculate mean
		fdJenk_m(i) = mean(fdJenk{i});

		% ------------------------------------------------------------------------------
		% Compute fd Power
		% ------------------------------------------------------------------------------
		% Threshold for flagging problem volumes
		fdPowerThr = 0.2;

		% Get FD
		fdPower{i} = GetFDPower(mov{i});

		% Calculate mean
		fdPower_m(i) = mean(fdPower{i});

		% ------------------------------------------------------------------------------
		% Get DVARS
		% ------------------------------------------------------------------------------
	    % Get dvars
	    dvars{i} = dlmread([workdir,'dvars.txt']);

	    % ------------------------------------------------------------------------------
	    % 1) Initial, gross movement exclusion
	    % ------------------------------------------------------------------------------
		% Calculate whether subject has suprathreshold mean movement
		% If the mean of displacement is greater than 0.55 mm (Sattethwaite), then exclude
		if fdJenk_m(i) > 0.55
			exclude(i,1) = 1;
		else
			exclude(i,1) = 0;
		end	

		% ------------------------------------------------------------------------------
		% 2) Stringent, multi criteria exclusion
		% ------------------------------------------------------------------------------
		% 2.1) Exclude on mean rms displacement
			% Calculate whether subject has suprathreshold mean movement
			% If the mean of displacement is greater than 0.2 mm (Ciric), then exclude
			if fdJenk_m(i) > 0.2
				x = 1;
			else
				x = 0;
			end	

		% 2.2) Exclude on proportion of spikes
			% Calculate whether subject has >20% suprathreshold spikes
			fdJenkThrPerc = round(numVols * 0.20);
			% If the number of volumes that exceed fdJenkThr are greater than %20, then exclude
			if sum(fdJenk{i} > fdJenkThr) > fdJenkThrPerc
				y = 1;
			else
				y = 0;
			end

		% 2.3) Exclude on large spikes (>5mm)
			if any(fdJenk{i} > 5)
				z = 1;
			else
				z = 0;
			end

		% If any of the above criteria are true of subject i, mark for exclusion
		if x == 1 | y == 1 | z == 1
			exclude(i,2) = 1;
		else
			exclude(i,2) = 0;
		end

	    % ------------------------------------------------------------------------------
	    % 3) Spike Regression
	    % ------------------------------------------------------------------------------
		% threshold for exclusion in minutes
		thresh = 4;

        % generate spike regressors
        spikereg = GetSpikeRegressors(fdJenk{i},fdJenkThr);

		% number of volumes - number of spike regressors (columns)
		numCVols = numVols - size(spikereg,2);

		% Compute length, in minutes, of time series data left after censoring
		NTime = (numCVols * TR)/60;

		% if less than threshold, mark for exclusion
		if NTime < thresh
			exclude(i,3) = 1;
		else
			exclude(i,3) = 0;
		end

	    % ------------------------------------------------------------------------------
	    % 4) JP12 Scrubbing
	    % ------------------------------------------------------------------------------
	    dvarsThr = 30;
		% scrubProximal = 'yes';
		scrubProximal = 'no';
		[JP12ScrubMask{i}, exclude(i,4)] = JP12_GetScrubMask(fdPower{i},dvars{i},TR,fdPowerThr,dvarsThr,scrubProximal);

	    % ------------------------------------------------------------------------------
	    % 5) JP14 Scrubbing
	    % ------------------------------------------------------------------------------
	    dvarsThr = 20;
		[scrubmask_temp, exclude(i,5)] = JP14_GetScrubMask(fdPower{i},dvars{i},TR,fdPowerThr,dvarsThr);
		JP14ScrubMask{i} = scrubmask_temp(:,2);	    

	    n = numel(msg);
	    fprintf(repmat('\b',1,n));
	end

	exclude = logical(exclude);
end
