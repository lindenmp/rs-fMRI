% ------------------------------------------------------------------------------
% This script simply loops over subjects and runs ThePlot
% 
% Copyright (C) 2017, Linden Parkes <lindenparkes@gmail.com>,
% ------------------------------------------------------------------------------
clear all; close all; clc

% ------------------------------------------------------------------------------
% Add paths - edit this section
% ------------------------------------------------------------------------------
funcdir = '/home/lindenmp/kg98/Linden/Scripts/rs-fMRI/func/';
addpath(funcdir)

% ------------------------------------------------------------------------------
% Set string switches
% ------------------------------------------------------------------------------
Projects = {'Beijing_Zang','UCLA','OCDPG','NYU_2'};
WhichProject = Projects{2};

noiseOptions = {'24P+8P',...
				'24P+8P+4GSR',...
				'24P+aCC',...
				'24P+aCC+4GSR',...
				'ICA-AROMA+2P',...
				'ICA-AROMA+2P+GSR',...
				'24P+8P+4GSR+SpikeReg',...
				'24P+4P+2GSR+JP14Scrub'};

noiseOptionsNames = {'24HMP+8Phys',...
					'24HMP+8Phys+4GSR',...
					'24HMP+aCompCor',...
					'24HMP+aCompCor+4GSR',...
					'ICA-AROMA+2Phys',...
					'ICA-AROMA+2Phys+GSR',...
					'24HMP+8Phys+4GSR+SpikeReg',...
					'24HMP+4Phys+2GSR+JP14Scrub'};

numPrePro = length(noiseOptions);

% ------------------------------------------------------------------------------
% Set project variables
% ------------------------------------------------------------------------------
% parentdir = '~/Dropbox/Work/ResProjects/';
parentdir = '/home/lindenmp/kg98_scratch/Linden/ResProjects/';
switch WhichProject
	case 'OCDPG'
		projdir = [parentdir,'rfMRI_denoise/OCDPG/'];
		sublist = [projdir,'M3_OCDPGe.csv'];
		datadir = [projdir,'data/'];
        preprostr = '/func/prepro/';
        t1str = '/anat/';

		TR = 2.5;
	case 'UCLA'
		projdir = [parentdir,'rfMRI_denoise/UCLA/'];
		sublist = [projdir,'UCLA.csv'];
		datadir = [projdir,'data/'];
        preprostr = '/func/prepro/';
        t1str = '/anat/';

		TR = 2;
	case 'NYU_2'
		projdir = [parentdir,'rfMRI_denoise/NYU_2/'];
		sublist = [projdir,'NYU_2.csv'];
		datadir = [projdir,'data/'];
		% Baseline data directory string
		% Note, we use the baseline data to calculate motion
		preprostr = '/session_1/func_1/prepro/';
		% preprostr = '/session_1/func_2/prepro/';
		% preprostr = '/session_2/func_1/prepro/';
        t1str = '/session_1/anat_1/';
	
		TR = 2;
	case 'COBRE'
		projdir = [parentdir,'rfMRI_denoise/COBRE/'];
		sublist = [projdir,'COBRE.csv'];
		datadir = [projdir,'data/'];
        preprostr = '/func/prepro/';
        t1str = '/anat/';

		TR = 2;
	case 'Beijing_Zang'
		projdir = [parentdir,'rfMRI_denoise/Beijing_Zang/'];
		sublist = [projdir,'Beijing_Zang.csv'];
		datadir = [projdir,'data/'];
        preprostr = '/func/prepro/';
        t1str = '/anat/';

		TR = 2;
end

% ------------------------------------------------------------------------------
% Setup output dir
% ------------------------------------------------------------------------------
% Setup output directory
outdir = ([projdir,'QC_ThePlot_Jenk']); 
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
% Subject list
% ------------------------------------------------------------------------------
fprintf(1, 'Loading metadata...\n');
metadata = readtable(sublist);

% ------------------------------------------------------------------------------
% Exclusion and censoring stuff
% ------------------------------------------------------------------------------
[metadata.exclude,metadata.mov,metadata.fdJenk,metadata.fdJenk_m,metadata.fdPower,metadata.fdPower_m,metadata.dvars,metadata.JP12ScrubMask,metadata.JP14ScrubMask] = GetExcludeForSample(datadir,metadata.ParticipantID,TR,preprostr);
fprintf(1, 'done\n');

% compute number of volumes using the length of fdJenk
% note, this is assumed to be same for all subjects!
numVols = length(metadata.fdJenk{1});

% compute numsubs
numSubs = size(metadata,1);

% ------------------------------------------------------------------------------
% Sort descending by mean rms displacement per subject
% ------------------------------------------------------------------------------
[srt,idx] = sort(metadata.fdJenk_m,'descend');
% [srt,idx] = sort(metadata.fdPower_m,'descend');

% ------------------------------------------------------------------------------
% Loop over noise correction methods
% ------------------------------------------------------------------------------
for n = 1:numPrePro
	WhichNoise = noiseOptions{n};
	WhichNoiseName = noiseOptionsNames{n};
	fprintf(1, '\nProcessing data: %s\n',WhichNoise);

	% ------------------------------------------------------------------------------
	% Loop over subjects in order of descending movement issues
	% ------------------------------------------------------------------------------
	% for i = 1:numSubs
	for i = [1,numSubs] % this will just do the highest and lowest motion subject
		fprintf(1,'\tProcessing subject %u/%u: %s\n',i,numSubs,metadata.ParticipantID{idx(i)})

		preprodir = [datadir,metadata.ParticipantID{idx(i)},preprostr,WhichNoise,'/'];
		t1dir = [datadir,metadata.ParticipantID{idx(i)},t1str];

		EPI = 'epi_prepro.nii';
		t1name = [metadata.ParticipantID{idx(i)},'_T1w.nii'];
		gmMask = ['gm50_bin.nii'];
		wmMask = ['wbc2c',t1name(1:end-4),'_e1.nii'];

	    filesToCheck = {{EPI}, ...
	                    {gmMask,wmMask}};

	    dirsOfFiles = {preprodir, ...
	                    t1dir};

        % Check if the epi_prepro.nii exists at all, if it doesnt, then the subject was processed
        % (this will happen if they were excluded according to optimized scrubbing procedures)
        if exist([dirsOfFiles{1},filesToCheck{1}{1},'.gz']) == 2 | exist([dirsOfFiles{1},filesToCheck{1}{1}]) == 2

		    for j = 1:length(filesToCheck)
		        for k = 1:length(filesToCheck{j})
		            if exist([dirsOfFiles{j},filesToCheck{j}{k},'.gz']) == 2
		                fprintf(1, '\t\tDecompressing %s\n',filesToCheck{j}{k});
		                gunzip([dirsOfFiles{j},filesToCheck{j}{k},'.gz'],dirsOfFiles{j})
		                delete([dirsOfFiles{j},filesToCheck{j}{k},'.gz'])
		            elseif exist([dirsOfFiles{j},filesToCheck{j}{k},'.gz']) == 0
		                fprintf(1, '\t\tNo need to decompress %s\n',filesToCheck{j}{k});
		            end
		        end
		    end

		    fprintf(1, '\t\tDrawing figure...\n');
			% ------------------------------------------------------------------------------
			% GetTSCompartment
			% ------------------------------------------------------------------------------
		    fsldir = '/usr/local/fsl/5.0.9/fsl/bin/'; % M3
			% fsldir = '/usr/share/fsl/5.0/bin/'; % local macbook
			[ts_compartment,key_compartment] = GetTSCompartment(fsldir,[preprodir,EPI],[t1dir,gmMask],[t1dir,wmMask]);
			% normalise
			ts_compartment = BF_NormalizeMatrix(ts_compartment,'maxmin');

			% ------------------------------------------------------------------------------
			% ThePlot
			% ------------------------------------------------------------------------------
			% Threshold for flagging problem volumes
			ThePlot([metadata.ParticipantID{idx(i)},' / ',WhichNoiseName],metadata.mov{idx(i)},metadata.fdPower{idx(i)},metadata.JP14ScrubMask{idx(i)},metadata.fdJenk{idx(i)},metadata.dvars{idx(i)},ts_compartment,key_compartment,TR)
			
			fig = gcf;
			set(fig,'PaperPositionMode','Auto')
			print(fig,[num2str(i),'_',metadata.ParticipantID{idx(i)},'_',WhichNoise,'.bmp'],'-dbmp')
			close all

			% Recompress files
		    for j = 1:length(filesToCheck)
		        for k = 1:length(filesToCheck{j})
		            fprintf(1, '\t\tRecompressing %s\n',filesToCheck{j}{k});
		            gzip([dirsOfFiles{j},filesToCheck{j}{k}],dirsOfFiles{j})
		            delete([dirsOfFiles{j},filesToCheck{j}{k}])
		        end
		    end
	    else
	    	fprintf(1, '\t\tData not present for %s / %s. Skipping...\n', metadata.ParticipantID{idx(i)}, WhichNoise);
	    end
	end
end
