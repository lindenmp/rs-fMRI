%% GetExcludeForSample: 
function [exclude,mov,fdJenk,fdJenk_m] = GetExcludeForSample(datadir,ParticipantIDs,str,mname)
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

	if nargin < 3
		str = '';
	end

	if nargin < 4
		mname = 'rp*txt';
	end

	numSubs = length(ParticipantIDs);

	% ------------------------------------------------------------------------------
	% Containers
	% ------------------------------------------------------------------------------
	mov = cell(numSubs,1);

	fdJenk = cell(numSubs,1);
	fdJenk_m = zeros(numSubs,1);

	exclude = zeros(numSubs,2);

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
	    mfile = dir([workdir,mname]);
	    mov{i} = dlmread([workdir,mfile(1).name]);

	    % ------------------------------------------------------------------------------
		% Compute fd (Jenkinson2002) (used by Satterthwaite2012)
	    % ------------------------------------------------------------------------------
		% Get FD
		fdJenk{i} = GetFDJenk(mov{i});

		% Calculate mean
		fdJenk_m(i) = mean(fdJenk{i});

	    % ------------------------------------------------------------------------------
	    % Initial, gross movement exclusion
	    % ------------------------------------------------------------------------------
		% 1) Exclude on mean rms displacement
 			% Calculate whether subject has suprathreshold mean movement
			% If the mean of displacement is greater than 0.55 mm (Sattethwaite), then exclude
			if fdJenk_m(i) > 0.55
				x = 1;
			else
				x = 0;
			end	

		% If any of the above criteria are true of subject i, mark for exclusion
		if x == 1
			exclude(i,1) = 1;
		else
			exclude(i,1) = 0;
		end

		clear x
		% ------------------------------------------------------------------------------
		% Stringent, multi criteria exclusion
		% ------------------------------------------------------------------------------
		% Threshold for detecting 'spikes'
		fdJenkThr = 0.25;

		% 1) Exclude on mean rms displacement
			% Calculate whether subject has suprathreshold mean movement
			% If the mean of displacement is greater than 0.2 mm (Ciric), then exclude
			if fdJenk_m(i) > 0.2
				x = 1;
			else
				x = 0;
			end	

		% 2) Exclude on proportion of spikes
			% Calculate whether subject has >20% suprathreshold spikes
			numVols = size(mov,1)-1;
			fdJenkThrPerc = round(numVols * 0.20);
			% If the number of volumes that exceed fdJenkThr are greater than %20, then exclude
			if sum(fdJenk{i} > fdJenkThr) > fdJenkThrPerc
				y = 1;
			else
				y = 0;
			end

		% 3) Exclude on large spikes (>5mm)
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

	    n = numel(msg);
	    fprintf(repmat('\b',1,n));
	end

	exclude = logical(exclude);
end