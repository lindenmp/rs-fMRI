function [fd] = GetFDVanD(mov)

	% This function compute framewise displacement according to Van Dijk's work
	%
	% ------
	% INPUTS
	% ------
	% mov       - an N x 3 matrix containing 3 translation movement parameters and where N
	%           = length of the time series.
	%
	% -------
	% OUTPUTS
	% -------
	% fd        - an N-1 length vector representing the total framewise
	%           displacement
	%
	% Linden Parkes, Brain & Mental Health Laboratory, 2016
	% ------------------------------------------------------------------------------
    
    fprintf(1,'Computing Van Dijk''s framewise displacement...');

	% detrend motion regressors
    mov = detrend(mov,'linear');

    % number of time points
    N = size(mov,1);

    % initialise fd variable.
    fd = zeros(N,1);
    
    % start at volume 2
    for i = 2:N
    	x = mov(i,1);
    	y = mov(i,2);
    	z = mov(i,3);

    	x_1 = mov(i-1,1);
    	y_1 = mov(i-1,2);
    	z_1 = mov(i-1,3);

    	fd(i) = rms([x y z]) - rms([x_1 y_1 z_1]);
    end

    fprintf(1,' done\n');

end