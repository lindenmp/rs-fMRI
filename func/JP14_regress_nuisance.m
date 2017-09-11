% ------------------------------------------------------------------------------
% Originally authored by Jonathan Power (see Power et al., 2014. NeuroImage)
% This code was downloaded from: http://www.jonathanpower.net/2014-ni-motion-2.html
% 
% Adapted by Linden Parkes for use in Parkes et al., 2017. bioRxiv
% Notes:
% 	this functions assumes that time is on the second dimension.
% 	But it assumes time is the first dimension on the noiseTS matrix.
% 	it also assumes that 1 = volumes to keep in tmask, not to discard.
% ------------------------------------------------------------------------------
function [data_out zb newregs] = JP14_regress_nuisance(data,noiseTS,tmask)

	[vox ts] = size(data);
	
	if nargin < 3;
		tmask = ones(ts,1);
	end

	% make logical
	tmask = logical(tmask);

	% ------------------------------------------------------------------------------
	% Regressors
	% ------------------------------------------------------------------------------
	% First, create regressors using the censored data
	zlinreg = noiseTS(tmask,:); % only censored data
	[zlinreg DMDTB] = JP14_demean_detrend(zlinreg'); % obtain fits for censored data
	zlinreg = zlinreg';
	zlinreg = zscore(zlinreg);

	% Next, use the fit from the above censored regressors to generate regressors that span the uncensored data too
	linreg = [repmat(1,[ts 1]) linspace(0,1,ts)'];
	newregs = DMDTB * linreg'; % predicted all regressors demean/detrend
	newregs = noiseTS - newregs'; % these are the demeaned detrended regressors
	newregs = zscore(newregs);

	% ------------------------------------------------------------------------------
	% fMRI data
	% ------------------------------------------------------------------------------
	% demean and detrend the censored data
	data_c = data(:,tmask);
	[data_c data_c_beta] = JP14_demean_detrend(data_c);

	% calculate betas on the censored data
	tempboldcell = num2cell(data_c',1);
	zlinregcell = repmat({zlinreg},[1 vox]);
	zb = cellfun(@mldivide,zlinregcell,tempboldcell,'uniformoutput',0);
	zb = cell2mat(zb);

	% demean and detrend all data using censored fits
	zmdttotimg = data_c_beta * linreg';
	zmdttotimg = data - zmdttotimg;

	% calculate residuals on all the data
	zb = zb';
	tempintvals = zb * newregs';
	data_out = zmdttotimg - tempintvals;
end
