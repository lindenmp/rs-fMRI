% ------------------------------------------------------------------------------
% Originally authored by Jonathan Power (see Power et al., 2014. NeuroImage).
% This code was downloaded from: http://www.jonathanpower.net/2014-ni-motion-2.html
% 
% Adapted by Linden Parkes for use in Parkes et al., 2017. bioRxiv
% This produces interpolated fMRI time series using censored data
% ------------------------------------------------------------------------------
function [H,f,s,c,tau,w] = JP14_getTransform(t,data,TH,Tr,ofac,hifac)

	% Inputs
	% 
	% t 	- column vector listing the *censored* time points in seconds
	%		e.g., t = ([1:N]') * TR;
	% 		2
	% 		4
	% 		6
	% 		12
	% 		14
	% 		16
	% 		Note the missing time points between 6 and 12
	% 
	% data 	- N (num time points) * numVoxels matrix of *censored* fMRI time series data
	% 
	% TH 	- column vector listing the *uncensored* time points in seconds
	% 
	% Tr 	- TR of fMRI sequences in seconds. e.g., TR = 2
	% 
	% ofac 	- oversampling frequency (generally >=4 )
	% hifac - highest frequency allowed

	% Outputs
	% H 	- the output time series reconstructed from the basis functions.
	% 		The censored sections of this time series need to be substituted back into the
	% 		original time series. That is, this is NOT the original time series + interpolated
	% 		volumes. It is simply the reconstructed time series.


	N = size(data,1); %Number of time points
	T = max(t) - min(t); %Total time span

	%calculate sampling frequencies
	f = (1/(T*ofac):1/(T*ofac):hifac*N/(2*T)).';

	%angular frequencies and constant offsets
	w = 2*pi*f;
	tau = atan2(sum(sin(2*w*t.'),2),sum(cos(2*w*t.'),2))./(2*w);

	%spectral power sin and cosine terms
	cterm = cos(w*t.' - repmat(w.*tau,1,length(t)));
	sterm = sin(w*t.' - repmat(w.*tau,1,length(t)));

	num_tseries = size(data,2); %Number of time series

	D = data;
	D = reshape(D,1,N,num_tseries);

	% This calculation is done by separately for the numerator, denominator,
	% and the division

	Cmult = bsxfun(@times, cterm,D);
	numerator = sum(Cmult,2); %Modify the numerator to get the cw term

	denominator = sum(bsxfun(@power,cterm,2),2); %use Cterm in place of cterm to make it exactly the denominator in the original expression
	C_final_new = bsxfun(@rdivide,numerator,denominator);
	c = C_final_new;
	clear numerator denominator cterm Cmult C_final_new

	%% Repeat the above for Sine term
	Smult = bsxfun(@times, sterm,D);
	numerator = sum(Smult,2); %Modify the numerator to get the sw term
	denominator = sum(bsxfun(@power,sterm,2),2);
	S_final_new = bsxfun(@rdivide,numerator,denominator);
	s = S_final_new;
	clear numerator denominator sterm Smult S_final_new

	% The inverse function to re-construct the original time series
	Time = TH';
	T_rep = repmat(Time,[size(f,1),1,size(data,2)]);
	% T_rep = bsxfun(@minus,T_rep,tau);
	% s(200) = s(199);
	w = 2*pi*f;
	prod = bsxfun(@times,T_rep,w);
	sin_t = sin(prod);
	cos_t = cos(prod);
	sw_p = bsxfun(@times,sin_t,s);
	cw_p = bsxfun(@times,cos_t,c);
	S = sum(sw_p);
	C = sum(cw_p);
	H = C + S;
	H = reshape(H,size(Time,2),size(data,2));

	%Normalize the reconstructed spectrum, needed when ofac > 1
	Std_H = std(H);
	Std_h = std(data);
	norm_fac = Std_H./Std_h;
	H = bsxfun(@rdivide,H,norm_fac);
end
