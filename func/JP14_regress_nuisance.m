% ------------------------------------------------------------------------------
% Originally authored by Jonathan Power (see Power et al., 2014. NeuroImage)
% This code was downloaded from: http://www.jonathanpower.net/2014-ni-motion-2.html
% 
% Adapted by Linden Parkes for use in Parkes et al., 2017. bioRxiv
% Notes:
% 	this functions assumes that time is on the second dimension.
% 	But it assumes time is the first dimension on the noiseTS matrix.
% 	it also assumes that 1 = volumes to keep in scrubmask, not to discard.
% ------------------------------------------------------------------------------
function [tempimg zb newregs] = regress_nuisance(tempimg,totregs,tot_tmask)

	[vox ts]=size(tempimg);
	zlinreg=totregs(logical(tot_tmask),:); % only desired data
	[zlinreg DMDTB]=demean_detrend(zlinreg'); % obtain fits for desired data
	zlinreg=zlinreg';
	zstd=std(zlinreg); % calculate std
	zmean=mean(zlinreg);
	zlinreg=zlinreg-(repmat(zmean,[size(zlinreg,1) 1]))./(repmat(zstd,[size(zlinreg,1) 1])); % zscore

	linreg=[repmat(1,[ts 1]) linspace(0,1,ts)'];
	newregs=DMDTB*linreg'; % predicted all regressors demean/detrend
	newregs=totregs-newregs'; % these are the demeaned detrended regressors
	newregs=newregs-(repmat(zmean,[size(newregs,1) 1]))./(repmat(zstd,[size(newregs,1) 1])); % zscore

	% now we have z-scored, detrended good and all regressors.

	% demean and detrend the desired data
	zmdtimg=tempimg(:,logical(tot_tmask));
	[zmdtimg zmdtbetas]=demean_detrend(zmdtimg);

	% calculate betas on the good data
	tempboldcell=num2cell(zmdtimg',1);
	zlinregcell=repmat({zlinreg},[1 vox]);
	zb = cellfun(@mldivide,zlinregcell,tempboldcell,'uniformoutput',0);
	zb=cell2mat(zb);

	% demean and detrend all data using good fits
	[zmdttotimg]=zmdtbetas*linreg';
	zmdttotimg=tempimg-zmdttotimg;

	% calculate residuals on all the data
	zb=zb';
	tempintvals=zb*newregs';
	tempimg=zmdttotimg-tempintvals;

end
