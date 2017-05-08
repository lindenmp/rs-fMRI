function [fd,delta_mov] = GetFDPower(mov, head)

	% This function compute framewise displacement according to Power's work
	%
	% ------
	% INPUTS
	% ------
	% mov       - an N x 6 matrix containing 6 movement parameters and where N
	%           = length of the time series.
	% head      - head radius (in mm) to use when converting radians to mm. 
	% 			default = 50mm, as in Power et al.
	%
	% -------
	% OUTPUTS
	% -------
	% fd        - an N-1 length vector representing the total framewise
	%           displacement
	% delta_mov - N-1 x 6 matrix of differentiated movement parameters
	%
	% Linden Parkes, Brain & Mental Health Laboratory, 2016
	% Alex Fornito, Brain & Mental Health Laboratory, 2016
	% ------------------------------------------------------------------------------

    if nargin < 2
    	head = 50;
    end

	% detrend motion regressors
    mov = detrend(mov,'linear');

	% convert degrees to radians to mm
	mov(:,4:6) = head*pi/180*mov(:,4:6);

	% differentiate movement parameters
    delta_mov = [
                zeros(1,size(mov,2));
				diff(mov);
				];

	% compute total framewise displacement
	fd = sum(abs(delta_mov),2);

end