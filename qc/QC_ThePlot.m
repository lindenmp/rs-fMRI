clear all; clc;

% ------------------------------------------------------------------------------
% This script simply loops over subjects and runs ThePlot
% 
% Linden Parkes, Brain & Mental Health Laboratory, 2016
% ------------------------------------------------------------------------------

% ------------------------------------------------------------------------------
% Select project
% ------------------------------------------------------------------------------
WhichProject = 'OCDPG'

switch WhichProject
	case 'OCDPG'
		outdir = '/gpfs/M2Home/projects/Monash076/Linden/OCDPG/';
		sublist = '/gpfs/M2Home/projects/Monash076/Linden/Sublists/OCDPG.txt';
        datadir = '/gpfs/M2Home/projects/Monash076/Linden/OCDPG/data/';

		% outdir = '~/Dropbox/Work/ResProjects/2014/OCDPG/'; 
		% sublist = '~/Dropbox/scripts/projects/OCDPG/sublists/SubjectIDs.txt';
		% datadir = '/media/lindenmp/ResProjects1/OCDPG/data/';
	case 'COBRE'
		outdir = '~/Dropbox/Work/ResProjects/2016/BF_TimeSeries/COBRE/'; 
		sublist = '/media/lindenmp/ResProjects2/COBRE/SubjectIDs.txt';
		datadir = '/media/lindenmp/ResProjects2/COBRE/data/';
	case 'UCLA'
		outdir = '~/Dropbox/Work/ResProjects/2016/BF_TimeSeries/UCLA/'; 
		sublist = '/media/lindenmp/ResProjects2/UCLA/SubjectIDs.txt';
		datadir = '/media/lindenmp/ResProjects2/UCLA/data/';
	case 'NAMIC'
		outdir = '~/Dropbox/Work/ResProjects/2016/BF_TimeSeries/NAMIC/'; 
		sublist = '/media/lindenmp/ResProjects2/NAMIC/SubjectIDs.txt';
		datadir = '/media/lindenmp/ResProjects2/NAMIC/data/';
end

% ------------------------------------------------------------------------------
% Load in Exclude data from qc_Exclude
% ------------------------------------------------------------------------------
cd(outdir)
load([WhichProject,'_Exclude.mat'])

numSubs = length(DarisIDs);

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
	subject = DarisIDs{idx(i)};
	fprintf(1,'Processing subject %u/%u: %s\n',i,numSubs,subject)

	switch WhichProject
		case 'OCDPG'
			rawdir = [datadir,subject,'/rfMRI/'];
			preprodir = [rawdir,'prepro/'];
			cleandir = [preprodir,'24P+aCC/'];
			
			t1dir = [datadir,subject,'/t1/'];
			
			rsData = 'srest_prepro.nii';
			dvarsExtract = 'idbwratepi.nii';
			
			movParams = 'rp_atepi.txt';

			t1name = 'ct1.nii';
			gmMask = ['wc1',t1name];
			wmMask = ['wc2',t1name];
			csfMask = ['wc3',t1name];
		
		case 'COBRE'
			rawdir = [datadir,subject,'/session_1/rest_1_Comp_GSR/'];
			preprodir = rawdir;
			cleandir = rawdir;

			t1dir = [datadir,subject,'/session_1/anat_1/'];

			rsData = 'srest_prepro.nii';
			dvarsExtract = 'detrend_4DVolume.nii';
			
			movParams = 'rp_arest.txt';

			t1name = 'mprage.nii';
			gmMask = ['crwc1',t1name];
			wmMask = ['crwc2',t1name];
			csfMask = ['crwc3',t1name];
		
		case 'UCLA'
			rawdir = [datadir,subject,'/func/'];
			preprodir = rawdir;
			cleandir = rawdir;

			t1dir = [datadir,subject,'/anat/'];
			
			rsData = 'srest_prepro.nii';
			
			movParams = ['rp_a',subject,'_task-rest_bold.txt'];

			t1name = [subject,'_T1w.nii'];
			gmMask = ['crwc1',t1name];
			wmMask = ['crwc2',t1name];
			csfMask = ['crwc3',t1name];
		
		case 'NAMIC'
			rawdir = [datadir,subject,'/rest_1/'];
			preprodir = rawdir;
			cleandir = rawdir;

			t1dir = [datadir,subject,'/anat_1/'];
			
			rsData = 'srest_prepro.nii';
			
			movParams = 'rp_arest.txt';
			
			t1name = 'mprage.nii';
			gmMask = ['crwc1',t1name];
			wmMask = ['crwc2',t1name];
			csfMask = ['crwc3',t1name];

	end

	% ------------------------------------------------------------------------------
	% Compute fd (Power2012)
	% ------------------------------------------------------------------------------
	% Threshold for flagging problem volumes
	fdPowerThr = 0.2;

	% Get FD
	fdPower{idx(i)} = GetFDPower(mov{idx(i)});

	% ------------------------------------------------------------------------------
	% Compute DVARS (Power2012)
	% ------------------------------------------------------------------------------
	dvars{idx(i)} = GetDVARS([preprodir,dvarsExtract],[t1dir,gmMask]);
	
	% ------------------------------------------------------------------------------
	% GetTSCompartment
	% ------------------------------------------------------------------------------
    fsldir = '/usr/local/fsl/5.0.9/bin/'; % directory where fsl is
	% fsldir = '/usr/share/fsl/5.0/bin/'; % directory where fsl is
	
	[ts_compartment,key_compartment] = GetTSCompartment(fsldir,[cleandir,rsData],[t1dir,gmMask],[t1dir,wmMask],[t1dir,csfMask]);

	% ------------------------------------------------------------------------------
	% ThePlot
	% ------------------------------------------------------------------------------
	ts_compartment = BF_NormalizeMatrix(ts_compartment,'maxmin');
	ThePlot(subject,mov{idx(i)},fdPower{idx(i)},fdJenk{idx(i)},dvars{idx(i)}/10,ts_compartment,key_compartment)

	cd(outdir)
	fig = gcf;
	pos = get(fig,'Position');
	
	set(fig,'PaperType','A4','PaperPositionMode','Auto','PaperUnits','CENTIMETERS','PaperPosition',[.63, .63, 19.72, 28.41])
	
	print(fig,[WhichProject,'ThePlot'],'-dpsc2','-append')
	close all
end