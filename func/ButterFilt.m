function [] = ButterFilt(FiltIn,MaskIn,TR,HighPass,LowPass)
	% This function performs 4th order Butterworth bandpass filtering on EPI ts
	% Note, requires read.m function from func/ folder
	%
	% ------
	% INPUTS
	% ------
	% FiltIn 	- location and name of processed resting state ts. e.g., /path/to/dir/ts.nii
	%
	% TR 		- TR of EPI data in seconds. e.g., 2
	% 
	% HighPass 	- Cut off for bottom end of bandpass in Hz. e.g., HighPass = 0.008;
	% 
	% LowPass 	- Cut off for top end of bandpass in Hz. e.g., LowPass = 0.08;
	% 
	% -------
	% OUTPUTS
	% -------
	% Filtered 4D volume named 'epi_filtered.nii'
	%
	% Linden Parkes, Brain & Mental Health Laboratory, 2016
	% ------------------------------------------------------------------------------

	fprintf(1, 'Running Butterworth filter...');
	% ------------------------------------------------------------------------------
	% Setup Butterworth filter
	% ------------------------------------------------------------------------------
	% Frequency
	Fs = 1/TR;
	% Nyquist
	Nq = Fs/2;
	% Bandpass cutoffs
	Wn = [HighPass LowPass]/Nq;
	
	% Build filter
	FiltOrder = 2;
	[b,a] = butter(FiltOrder,Wn,'bandpass');

	% ------------------------------------------------------------------------------
	% Load in EPI ts
	% ------------------------------------------------------------------------------
	% Load in ts and a single time series
	[hdr,ts] = read(FiltIn);

	% Convert to double just incase ts comes in as something else
	% butter/filtfilt in matlab don't work on anything other than double ts
	ts = double(ts);

	% Find dimensions of 4D EPI matrix
	dim = size(ts);
	N = dim(4);

	% reshape to 2d
	ts = reshape(ts,[],N);
	% Put time on first dimension
	ts = double(ts');

	% ------------------------------------------------------------------------------
	% Load in brain mask
	% ------------------------------------------------------------------------------
	[hdr_mask,mask] = read(MaskIn);

	% reshape
	mask = reshape(mask,[],1);
	mask = mask'; mask = logical(mask);

	maskIdx = find(mask > 0);

	% ------------------------------------------------------------------------------
	% Demean
	% ------------------------------------------------------------------------------
	m = mean(ts,1);
	for i = 1:N
		ts(i,:) = ts(i,:) - m;
	end

	% ------------------------------------------------------------------------------
	% Run filter
	% ------------------------------------------------------------------------------
	ts_out = zeros(size(ts));

	ts_out(:,maskIdx) = filtfilt(b,a,ts(:,maskIdx));

	% ------------------------------------------------------------------------------
	% Add mean back
	% ------------------------------------------------------------------------------
	for i = 1:N
		ts_out(i,:) = ts_out(i,:) + m;
	end	

	% ------------------------------------------------------------------------------
	% Write out filtered ts
	% ------------------------------------------------------------------------------
	% Put time back on second dimensions
	ts_out = ts_out';
	% Reshape back to 4D matrix
	ts_out = reshape(ts_out,dim);

	% Write
	write(hdr,ts_out,'epi_filtered.nii')

	fprintf(1, 'done\n');

end