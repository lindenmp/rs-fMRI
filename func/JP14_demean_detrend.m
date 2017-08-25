% ------------------------------------------------------------------------------
% Originally authored by Jonathan Power (see Power et al., 2014. NeuroImage).
% This code was downloaded from: http://www.jonathanpower.net/2014-ni-motion-2.html
% 
% Adapted by Linden Parkes for use in Parkes et al., 2017. bioRxiv
% This produces identical output to matlab's detrend function except also outputs the betas
% ------------------------------------------------------------------------------
function [tempbold tempbetas] = demean_detrend(img,varargin)

	if ~isnumeric(img)
	    [tempbold]=read_4dfpimg(img); % read image
	else
	    [tempbold]=img;
	    clear img;
	end
	[vox ts]=size(tempbold);

	if ~isempty(varargin)
	    tmask=varargin{1,1};
	else
	    tmask=ones(ts,1);
	end

	linreg=[repmat(1,[ts 1]) linspace(0,1,ts)'];
	tempboldcell=num2cell(tempbold(:,logical(tmask))',1);
	linregcell=repmat({linreg(logical(tmask),:)},[1 vox]);
	tempbetas = cellfun(@mldivide,linregcell,tempboldcell,'uniformoutput',0);
	tempbetas=cell2mat(tempbetas);
	tempbetas=tempbetas';
	tempintvals=tempbetas*linreg';
	tempbold=tempbold-tempintvals;

	if nargin==3
	    outname=varargin{1,2};
	    write_4dfpimg(tempbold,outname,'bigendian');
	    write_4dfpifh(outname,size(tempbold,2),'bigendian');
	end
end
