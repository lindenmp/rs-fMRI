% ------------------------------------------------------------------------------
% This script simply loops over subjects and runs ThePlot
% 
% Linden Parkes, Brain & Mental Health Laboratory, 2016
% ------------------------------------------------------------------------------
clear all; close all; clc

% ------------------------------------------------------------------------------
% Add paths - edit this section
% ------------------------------------------------------------------------------
WhichMASSIVE = 'M3';
switch WhichMASSIVE
    case 'M2'
		funcdir = '/gpfs/M2Home/projects/Monash076/Linden/scripts/rs-fMRI/func/';
		addpath(funcdir)
		funcdir1 = '/gpfs/M2Home/projects/Monash076/Linden/scripts/Func/';
		addpath(funcdir1)
    case 'M3'
		funcdir = '/home/lindenmp/kg98/Linden/Scripts/rs-fMRI/func/';
		addpath(funcdir)
		funcdir1 = '/home/lindenmp/kg98/Linden/Scripts/Func/';
		addpath(funcdir1)
end

% ------------------------------------------------------------------------------
% Set string switches
% ------------------------------------------------------------------------------
Projects = {'OCDPG','UCLA','NYU_2','M3_COBRE','M3_UCLA','M3_NAMIC'};
WhichProject = Projects{4};

WhichNoise = 'sICA-AROMA+2P';

% ------------------------------------------------------------------------------
% Set project variables
% ------------------------------------------------------------------------------
switch WhichProject
	case 'OCDPG'
		projdir = '/gpfs/M2Home/projects/Monash076/Linden/OCDPG/';
		sublist = [projdir,'OCDPG.txt'];
		datadir = [projdir,'data/'];
		preprostr = '/rfMRI/prepro/';
		t1str = '/t1/';

		EPI = 'epi_prepro.nii';
		EPIdvars = 'idbwrdatepi.nii';
		movParams = 'rp_datepi.txt';
		t1name = 'ct1.nii';
		gmMask = ['wc1',t1name];
		wmMask = ['wc2',t1name];
		csfMask = ['wc3',t1name];

		TR = 2.5;
	case 'UCLA'
		projdir = '/gpfs/M2Home/projects/Monash076/Linden/UCLA/';
		sublist = [projdir,'UCLA.txt'];
		datadir = [projdir,'data/'];
        preprostr = '/rfMRI/prepro/';
		t1str = '/t1/';

		EPI = 'epi_prepro.nii';
		EPIdvars = 'idbwrdatepi.nii';
		movParams = 'rp_datepi.txt';
		t1name = 'ct1.nii';
		gmMask = ['wc1',t1name];
		wmMask = ['wc2',t1name];
		csfMask = ['wc3',t1name];

		TR = 2;
	case 'NYU_2'
		projdir = '/gpfs/M2Home/projects/Monash076/Linden/NYU_2/';
		sublist = [projdir,'NYU_2.txt'];
		datadir = [projdir,'data/'];
		preprostr = '/session_1/rest_1/prepro/';
		% preprostr = '/session_1/rest_2/prepro/';
		% preprostr = '/session_2/rest_1/prepro/';
		t1str = '/t1/';

		EPI = 'epi_prepro.nii';
		EPIdvars = 'idbwrdatrest.nii';
		movParams = 'rp_datrest.txt';
		t1name = 'anat.nii';
		gmMask = ['wc1',t1name];
		wmMask = ['wc2',t1name];
		csfMask = ['wc3',t1name];

		TR = 2;
	case 'M3_COBRE'
		projdir = '/home/lindenmp/kg98/Linden/ResProjects/SCZ_HCTSA/COBRE/';
		sublist = [projdir,'COBRE_SubjectIDs.txt'];
		datadir = [projdir,'data/'];
		preprostr = '/session_1/rest_1/prepro/';
		t1str = '/session_1/anat_1/';

		EPI = 'epi_prepro.nii';
		EPIdvars = 'idbwrdatrest.nii';
		movParams = 'rp_datrest.txt';
		t1name = 'cmprage.nii';
		gmMask = ['wc1',t1name];
		wmMask = ['wc2',t1name];
		csfMask = ['wc3',t1name];

		TR = 2;
	case 'M3_UCLA'
		projdir = '/home/lindenmp/kg98/Linden/ResProjects/SCZ_HCTSA/UCLA/';
		sublist = [projdir,'UCLA_SubjectIDs.txt'];
		datadir = [projdir,'data/'];
		preprostr = '/func/prepro/';
		t1str = '/anat/';

		TR = 2;
	case 'M3_NAMIC'
		projdir = '/home/lindenmp/kg98/Linden/ResProjects/SCZ_HCTSA/NAMIC/';
		sublist = [projdir,'NAMIC_SubjectIDs.txt'];
		datadir = [projdir,'data/'];
		preprostr = '/func/prepro/';
		t1str = '/anat/';

		TR = 3;
end

% ------------------------------------------------------------------------------
% Setup output dir
% ------------------------------------------------------------------------------
% Setup output directory
outdir = ([projdir,'QC_ThePlot']); 
if exist(outdir) == 0
	fprintf(1,'Initialising outdir\n')
	mkdir(outdir)
elseif exist(outdir) == 7
	fprintf(1,'Cleaning and re-initialising outdir\n')
	rmdir(outdir,'s')
	mkdir(outdir)
end

cd(outdir)

% ------------------------------------------------------------------------------
% Subject list
% ------------------------------------------------------------------------------
fileID = fopen(sublist);
ParticipantIDs = textscan(fileID,'%s');
ParticipantIDs = ParticipantIDs{1};

% compute numsubs
numSubs = length(ParticipantIDs);

% ------------------------------------------------------------------------------
% Jenkinson's mean FD
% ------------------------------------------------------------------------------
fprintf(1, 'Loading Jenkinson''s mean FD metric\n');
[~,mov,fdJenk,fdJenk_m] = GetExcludeForSample(datadir,ParticipantIDs,preprostr);
fprintf(1, 'done\n');

% compute number of volumes using the length of fdJenk
% note, this is assumed to be same for all subjects!
numVols = length(fdJenk{1});

% ------------------------------------------------------------------------------
% Sort descending by mean rms displacement per subject
% ------------------------------------------------------------------------------
[srt,idx] = sort(fdJenk_m,'descend');

fdPower = cell(numSubs,1);
dvars = cell(numSubs,1);

% ------------------------------------------------------------------------------
% Loop over subjects in order of descending movement issues
% ------------------------------------------------------------------------------
for i = 1:numSubs
	fprintf(1,'Processing subject %u/%u: %s\n',i,numSubs,ParticipantIDs{idx(i)})

	preprodir = [datadir,ParticipantIDs{idx(i)},preprostr];
	cleandir = [preprodir,WhichNoise,'/'];
	t1dir = [datadir,ParticipantIDs{idx(i)},t1str];

	if ismember('M3_UCLA',WhichProject,'rows') | ismember('M3_NAMIC',WhichProject,'rows')
		EPI = 'epi_prepro.nii';
		EPIdvars = ['idbwrdat',ParticipantIDs{idx(i)},'_task-rest_bold.nii'];
		movParams = ['rp_dat',ParticipantIDs{idx(i)},'_task-rest_bold.txt'];
		t1name = [ParticipantIDs{idx(i)},'_T1w.nii'];
		gmMask = ['wc1',t1name];
		wmMask = ['wc2',t1name];
		csfMask = ['wc3',t1name];
	end

    % Decompress files
    filesToCheck = {EPI, ...
                    EPIdvars, ...
                    gmMask, ...
                    wmMask, ...
                    csfMask};

    dirsOfFiles = {cleandir, ...
                    preprodir, ...
                    t1dir, ...
                    t1dir, ...
                    t1dir};

    for j = 1:length(filesToCheck)
        if exist([dirsOfFiles{j},filesToCheck{j},'.gz']) == 2
            fprintf(1, '\t\t Decompressing %s\n',filesToCheck{j});
            gunzip([dirsOfFiles{j},filesToCheck{j},'.gz'],dirsOfFiles{j})
            delete([dirsOfFiles{j},filesToCheck{j},'.gz'])
        elseif exist([dirsOfFiles{j},filesToCheck{j},'.gz']) == 0
            fprintf(1, '\t\t No need to decompress %s\n',filesToCheck{j});
        end
    end

	% Get fdPower
	fdPower{idx(i)} = GetFDPower(mov{idx(i)});

	% Compute DVARS
	dvars{idx(i)} = GetDVARS([preprodir,EPIdvars],[t1dir,gmMask]);
	
	% GetTSCompartment
    fsldir = '/usr/local/fsl/5.0.9/fsl/bin/'; % M3
    % fsldir = '/usr/local/fsl/5.0.9/bin/'; % M2
	% fsldir = '/usr/share/fsl/5.0/bin/'; % local macbook
	[ts_compartment,key_compartment] = GetTSCompartment(fsldir,[cleandir,EPI],[t1dir,gmMask],[t1dir,wmMask],[t1dir,csfMask]);
	% normalise
	ts_compartment = BF_NormalizeMatrix(ts_compartment,'maxmin');

	% ------------------------------------------------------------------------------
	% ThePlot
	% ------------------------------------------------------------------------------
	% Threshold for flagging problem volumes
	ThePlot(ParticipantIDs{idx(i)},mov{idx(i)},fdPower{idx(i)},fdJenk{idx(i)},dvars{idx(i)}/10,ts_compartment,key_compartment,TR)
	
	fig = gcf;
	set(fig,'PaperPositionMode','Auto')
	print(fig,[num2str(i),'_',ParticipantIDs{idx(i)},'.bmp'],'-dbmp')
	close all

	% Recompress files
    for j = 1:length(filesToCheck)
        fprintf(1, '\t\t Recompressing %s\n',filesToCheck{j});
        gzip([dirsOfFiles{j},filesToCheck{j}],dirsOfFiles{j})
        delete([dirsOfFiles{j},filesToCheck{j}])
    end
end
