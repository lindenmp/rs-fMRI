function [fd] = GetFDJenk(mov, head)

	% This function compute framewise displacement according to Jenkinson
	%
	% ------
	% INPUTS
	% ------
	% mov       - an N x 6 matrix containing 6 movement parameters and where N
	%           = length of the time series.
	% head      - head radius (in mm). default = 80mm, as in Yan/FSL.
	%
	% -------
	% OUTPUTS
	% -------
	% fd        - an N-1 length vector representing the total framewise
	%           displacement
	%
	% Linden Parkes, Brain & Mental Health Laboratory, 2016
	% Ben Fulcher, Brain & Mental Health Laboratory, 2016
	% ------------------------------------------------------------------------------
    
    if nargin < 2
    	head = 80;
    end

    % detrend motion regressors
    mov = detrend(mov,'linear');
    
    % number of time points
    N = size(mov,1);

    % initialise fd variable.
    fd = zeros(N,1);

    % start at volume 2
    for i = 2:N

    	% rigid body transform mat from time point i
    	T = GetTMat(mov(i,:));

    	% rigid body transform mat from time point i-1
    	T_1 = GetTMat(mov(i-1,:));

    	Y = T * inv(T_1) - eye(size(T));

    	A = Y(1:3,1:3);
    	b = Y(1:3,4);

    	x = mov(i,1); y = mov(i,2); z = mov(i,3);
    	m = b + A * [x; y; z];

    	% store in i-1 slot. fd(1) = displacement difference betwen timepoint 2 and timepoint 1
    	fd(i) = sqrt(1/5 * head^2 * trace(A'*A) + m'*m);
    end

end