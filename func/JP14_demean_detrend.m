% ------------------------------------------------------------------------------
% Originally authored by Jonathan Power (see Power et al., 2014. NeuroImage).
% This code was downloaded from: http://www.jonathanpower.net/2014-ni-motion-2.html
% 
% Adapted by Linden Parkes for use in Parkes et al., 2017. bioRxiv
% This produces identical output to matlab's detrend function except also outputs the betas
% ------------------------------------------------------------------------------
function [tempbold,tempbetas] = JP14_demean_detrend(data,scrubmask);

    tempbold = data;
	[vox ts] = size(tempbold); % Note, Power assumes time on second dimension, whereas we work with time on first dimension.

	if nargin < 2
		% Note that in this scrubmask, 1 = volumes to keep, not volumes to censor!
	    scrubmask = ones(ts,1);
	end
	
	% make logical
	scrubmask = logical(scrubmask);

	linreg = [repmat(1,[ts 1]) linspace(0,1,ts)'];
	tempboldcell = num2cell(tempbold(:,scrubmask)',1);
	linregcell = repmat({linreg(scrubmask,:)},[1 vox]);
	tempbetas = cellfun(@mldivide,linregcell,tempboldcell,'uniformoutput',0);
	tempbetas = cell2mat(tempbetas);
	tempbetas = tempbetas';
	tempintvals = tempbetas * linreg';
	tempbold = tempbold - tempintvals;

end
