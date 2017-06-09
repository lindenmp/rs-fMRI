% ------------------------------------------------------------------------------
% This script compiles time series and metadata for the SCZ_HCTSA project
% 
% Linden Parkes, Brain & Mental Health Laboratory, 2017
% ------------------------------------------------------------------------------
clear all; close all; clc

% ------------------------------------------------------------------------------
% Set string switches
% ------------------------------------------------------------------------------
Projects = {'OCDPG','UCLA','NYU_2','M3_COBRE','M3_UCLA','M3_NAMIC'};
WhichProject = Projects{4};

removeNoise = 'sICA-AROMA+2P+GSR';

WhichParc = 'Gordon';

% ------------------------------------------------------------------------------
% Set project variables
% ------------------------------------------------------------------------------
switch WhichProject
	case 'OCDPG'
		projdir = '~/Dropbox/Work/ResProjects/rfMRI_denoise/OCDPG/';
		sublist = [projdir,'OCDPGe.csv'];
		datadir = [projdir,'data/'];
		preprostr = '/rfMRI/prepro/';

		TR = 2.5;
	case 'UCLA'
		projdir = '~/Dropbox/Work/ResProjects/rfMRI_denoise/UCLA/';
		sublist = [projdir,'UCLA.csv'];
		datadir = [projdir,'data/'];
        preprostr = '/rfMRI/prepro/';

		TR = 2;
	case 'NYU_2'
		projdir = '~/Dropbox/Work/ResProjects/rfMRI_denoise/NYU_2/';
		sublist = [projdir,'NYU_2.csv'];
		datadir = [projdir,'data/'];
		% Baseline data directory string
		% Note, we use the baseline data to calculate motion
		preprostr = '/session_1/rest_1/prepro/';
		% preprostr = '/session_1/rest_2/prepro/';
		% preprostr = '/session_2/rest_1/prepro/';
	
		TR = 2;
	case 'M3_COBRE'
		projdir = '~/Dropbox/Work/ResProjects/SCZ_HCTSA/COBRE/';
		sublist = [projdir,'COBRE.csv'];
		datadir = [projdir,'data/'];
		preprostr = '/session_1/rest_1/prepro/';
	
		TR = 2;
	case 'M3_UCLA'
		projdir = '~/Dropbox/Work/ResProjects/SCZ_HCTSA/UCLA/';
		sublist = [projdir,'UCLA.csv'];
		datadir = [projdir,'data/'];
		preprostr = '/func/prepro/';
	
		TR = 2;
	case 'M3_NAMIC'
		projdir = '~/Dropbox/Work/ResProjects/SCZ_HCTSA/NAMIC/';
		sublist = [projdir,'NAMIC.csv'];
		datadir = [projdir,'data/'];
		preprostr = '/func/prepro/';
	
		TR = 3;
end

% ------------------------------------------------------------------------------
% Set parcellation
% Note, code is not setup to process multiple parcellations concurrently.
% ------------------------------------------------------------------------------
switch WhichParc
	case 'Gordon'
		Parc = 1;
	case 'Power'
		Parc = 2;
end

% ------------------------------------------------------------------------------
% Non-imaging metadata
% ------------------------------------------------------------------------------
metadata = readtable(sublist);
% convert participant IDs to strings
if ~iscellstr(metadata.ParticipantID)
	metadata.ParticipantID =  cellfun(@num2str, num2cell(metadata.ParticipantID), 'UniformOutput', false);
	if ismember('M3_COBRE',WhichProject,'rows')
		metadata.ParticipantID = cellfun(@(c)['00' c],metadata.ParticipantID,'UniformOutput',false);
	end
end

numSubs = size(metadata,1);

roiTS = cell(numSubs,1);

% ------------------------------------------------------------------------------
% Load in time series data
% ------------------------------------------------------------------------------
for i = 1:numSubs
	workdir = [datadir,metadata.ParticipantID{i},preprostr,removeNoise,'/'];

	clear temp
	temp = load([workdir,'cfg.mat']);

	cfg = temp.cfg;

	roiTS{i} = cfg.roiTS{Parc};
end

% ------------------------------------------------------------------------------
% Jenkinson's mean FD
% ------------------------------------------------------------------------------
fprintf(1, 'Loading Jenkinson''s mean FD metric\n');
[exclude,~,fdJenk,fdJenk_m] = GetExcludeForSample(datadir,metadata.ParticipantID,preprostr);
fprintf(1, 'done\n');

% compute number of volumes using the length of fdJenk
% note, this is assumed to be same for all subjects!
numVols = length(fdJenk{1});

% ------------------------------------------------------------------------------
% Create exclusion vector based on <4 minutes of post-scrub BOLD data
% Note, subjects not excluded yet. that is done during the preprocessing loop below
% ------------------------------------------------------------------------------
% exclusion vector
excludeCensor = zeros(numSubs,1);
% threshold for exclusion in minutes
thresh = 4;

% Loop over subjects and check if censoring would reduce data to <4 minutes
for i = 1:numSubs
	% Get spike regressors for subject i
	spikereg = GetSpikeRegressors(fdJenk{i},0.25);
	% number of volumes - number of spike regressors (columns)
	numCVols = numVols - size(spikereg,2);

	% Compute length, in minutes, of time series data left after censoring
	NTime = (numCVols * TR)/60;
	% if less than threshold, mark for exclusion
	if NTime < thresh
		excludeCensor(i) = 1;
	end
end
% convert to logical
excludeCensor = logical(excludeCensor);

exclude = [exclude excludeCensor];

fprintf(1, '%u subjects marked for exclusion based on >.55mm mean FD \n', sum(exclude(:,1)));
fprintf(1, '%u subjects marked for exclusion based on multi FD criteria \n', sum(exclude(:,2)));
fprintf(1, '%u subjects marked for exclusion based on 4 minute censoring criteria \n', sum(exclude(:,3)));

% ------------------------------------------------------------------------------
% Save
% ------------------------------------------------------------------------------
cd(projdir)

save([WhichProject,'_',removeNoise,'_',WhichParc,'.mat'],'metadata','roiTS','exclude')
