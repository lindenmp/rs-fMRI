function [] = SpatialNormalisationANTs(EPI,N,meanEPI,t1file,gm,wm,csf,mni_template,antsdir,funcdir,fUpfile)
	% This function will perform the following spatial pre-processing steps using ANTs methods
	% 1 - coregistration of mean EPI to T1
	% 2 - spatial normalization of T1 to MNI
	% 3 - Application of the normalization parameters to coregistered EPI
	% NOTE, this script requires that RealignEPI.m has been run before hand
	% ------
	% INPUTS
	% ------
	% EPI 			- path and base filename (i.e., without .nii extension) of 3D EPI files
	%          		e.g., '/path/to/epi/ratepi'
	% 				This is done because ANTs does not seem to like 4D files - atleast, I couldnt get it to work on MASSIVE! :\
	% 
	% N 			- number of volumes in EPI data
	% 
	% meanEPI 		- path and filename of a mean 3D input EPI file output from realignment
	%           	e.g., '/path/to/epi/meanatepi.nii'
	% 
	% t1file 		- path and filename of a native t1 image
	%           	e.g., '/path/to/t1/t1.nii'
	% 
	% gm 			- path and filename of a grey matter probability mask in native space (e.g., from SPMs New Segment routine)
	% wm 			- path and filename of a white matter probability mask in native space (e.g., from SPMs New Segment routine)
	% csf 			- path and filename of a csf probability mask in native space (e.g., from SPMs New Segment routine)
	% 
	% mni_template 	- path and filename of a template used for spatial normalization
	% 				e.g., [spmdir,'templates/T1.nii'];
	% 
	% antsdir 		- path to ants binaries
	% 
	% fundir 		- path to rfMRI-Func dir 
	% 				this is only because ANTs includes some scripts which arent part of their binary package, so I store them in the rfMRI-Func dir on GitHub
	% 
	% fUpfile		- (optional) path to an addition T1 acquired face up which should be input if t1file was acquired face down.  
	%               
	% -------
	% OUTPUTS
	% -------
	%
	% 
	% Linden Parkes, Brain & Mental Health Laboratory, 2016
	% ------------------------------------------------------------------------------

	ANTsCall = 'antsRegistrationSyN.sh';
	% ANTsCall = 'antsRegistrationSyNQuick.sh'; % a quicker variant of ANTs (fewer iterations)

	ANTsCall2 = 'antsRegistrationSyNFDown.sh';
	% ANTsCall2 = 'antsRegistrationSyNFDownQuick.sh'; % a quicker variant of ANTs (fewer iterations)

	% ------------------------------------------------------------------------------------
	% Rigid body coreg of mean epi to native T1 image
	% Note, we only ant the .mat file from this step
	% ------------------------------------------------------------------------------------
	movImage = meanEPI;
	refImage = t1file;
	output = 'epi2t1';
	transform = 'r';

	system([funcdir,ANTsCall,' -d 3 -m ',movImage,' -f ',refImage,' -t ',transform,' -o ',output]);

	delete('epi2t1.nii.gz')

	if nargin == 11
		% ------------------------------------------------------------------------------
		% Rigid body coreg fDown to fUP
		% ------------------------------------------------------------------------------
		movImage = t1file;
		refImage = fUpfile;
		output = 'fDown2fUp';
		transform = 'r';

		system([funcdir,ANTsCall,' -d 3 -m ',movImage,' -f ',refImage,' -t ',transform,' -o ',output]);

		delete('fDown2fUp.nii.gz')
	end

	% ------------------------------------------------------------------------------
	% Nonlinear warp of native T1 image to SPM template (i.e., MNI space)
	% ------------------------------------------------------------------------------
	movImage = t1file;
	refImage = mni_template;
	output = 't12MNI';
	transform = 's';

	if nargin < 11
		system([funcdir,ANTsCall,' -d 3 -m ',movImage,' -f ',refImage,' -t ',transform,' -o ',output]);
	elseif nargin == 11
		system([funcdir,ANTsCall2,' -d 3 -m ',movImage,' -f ',refImage,' -t ',transform,' -o ',output]);
	end

	% unzip normalised T1
	gunzip([output,'.nii.gz'])
	delete([output,'.nii.gz'])

	% % rename output T1 file
	[fPath,fName,fExt] = fileparts(movImage);
	warpedT1 = ['w',fName,fExt];
	movefile('t12MNI.nii',warpedT1)

	% ------------------------------------------------------------------------------
	% Apply transforms
	% ------------------------------------------------------------------------------
	refImage = mni_template;

	% The ANTs documentation says that the transforms are concatenated in the order they are entered into the command line
	% But! when I do it the command line output indicates that they are applied in the reverse to how they are entered
	% So, I enter them in reverse here...
	% warps = '-t epi2t1_0GenericAffine.mat -t t12MNI_0GenericAffine.mat -t t12MNI_1Warp.nii.gz ';
	warps = '-t t12MNI_1Warp.nii.gz -t t12MNI_0GenericAffine.mat -t epi2t1_0GenericAffine.mat ';

	% 1) 4D EPI
	for i = 1:N

	    if i < 10
	        ZeroPad = '0000';
	    elseif i >= 10 & i < 100
	        ZeroPad = '000';
	    elseif i >= 100
	        ZeroPad = '00';
	    end

		movImage = [EPI,'_',ZeroPad,num2str(i),'.nii'];
		[fPath,fName,fExt] = fileparts(movImage);
		output = ['w',fName,fExt];

		system([antsdir,'antsApplyTransforms -d 3 -e 0 -i ',movImage,' -r ',refImage,' -o ',output,' -n Linear ',warps]);

	end

	% 2) mean EPI
	movImage = meanEPI;
	[fPath,fName,fExt] = fileparts(movImage);
	output = ['w',fName,fExt];

	system([antsdir,'antsApplyTransforms -d 3 -e 0 -i ',movImage,' -r ',refImage,' -o ',output,' -n Linear ',warps]);

	% redefine warps
	% warps = '-t t12MNI_0GenericAffine.mat -t t12MNI_1Warp.nii.gz ';
	warps = '-t t12MNI_1Warp.nii.gz -t t12MNI_0GenericAffine.mat ';

	% 3) GM
	movImage = gm;
	[fPath,fName,fExt] = fileparts(movImage);
	output = ['w',fName,fExt];

	system([antsdir,'antsApplyTransforms -d 3 -e 0 -i ',movImage,' -r ',refImage,' -o ',output,' -n Linear ',warps]);

	% 4) WM
	if iscell(wm)
		for i = 1:length(wm)
			movImage = wm{i};
			[fPath,fName,fExt] = fileparts(movImage);
			output = ['w',fName,fExt];

			system([antsdir,'antsApplyTransforms -d 3 -e 0 -i ',movImage,' -r ',refImage,' -o ',output,' -n Linear ',warps]);
		end
	elseif ~iscell(wm)
		movImage = wm;
		[fPath,fName,fExt] = fileparts(movImage);
		output = ['w',fName,fExt];

		system([antsdir,'antsApplyTransforms -d 3 -e 0 -i ',movImage,' -r ',refImage,' -o ',output,' -n Linear ',warps]);
	end

	% 5) CSF
	if iscell(csf)
		for i = 1:length(csf)
			movImage = csf{i};
			[fPath,fName,fExt] = fileparts(movImage);
			output = ['w',fName,fExt];

			system([antsdir,'antsApplyTransforms -d 3 -e 0 -i ',movImage,' -r ',refImage,' -o ',output,' -n Linear ',warps]);
		end
	elseif ~iscell(csf)
		movImage = csf;
		[fPath,fName,fExt] = fileparts(movImage);
		output = ['w',fName,fExt];

		system([antsdir,'antsApplyTransforms -d 3 -e 0 -i ',movImage,' -r ',refImage,' -o ',output,' -n Linear ',warps]);
	end

end
