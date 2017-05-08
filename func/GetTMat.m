%% GetTMat: 
function [T] = GetTMat(mov)
	% This function produces a rigid body transform matrix for a given timepoint from 
	% movement parameters derived from SPM8's realignment step
	%
	% ------
	% INPUTS
	% ------
	% mov       - an 1 x 6 vector containing 6 movement parameters for a SINGLE time point
	%
	% -------
	% OUTPUTS
	% -------
	% T         - rigid body transormation matrix, T
	%
	% Linden Parkes, Brain & Mental Health Laboratory, 2016
	% Ben Fulcher, Brain & Mental Health Laboratory, 2016
	% ------------------------------------------------------------------------------

	x = mov(1,1);
	y = mov(1,2);
	z = mov(1,3);

	a = mov(1,4);
	b = mov(1,5);
	g = mov(1,6);

	% translation matrix
	t1 = [eye(3), [x; y; z]; zeros(1,3), 1];

	% alpha rotation matrix
	t2 = [1 0 0 0;...
		0 cos(a) sin(a) 0;...
		0 -sin(a) cos(a) 0;...
		0 0 0 1];

	% beta rotation matrix
	t3 = [cos(b) 0 sin(b) 0;...
		0 1 0 0;...
		-sin(b) 0 cos(b) 0;...
		0 0 0 1];

	% gamma rotation matrix
	t4 = [cos(g) sin(g) 0 0;...
		-sin(g) cos(g) 0 0;...
		0 0 1 0;...
		0 0 0 1];

	T = t1*t2*t3*t4;
end