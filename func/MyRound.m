%% MyRound: the inbuilt round function on MASSIVE's copy of matlab doesnt work
% Gives a 'too many inputs' error.
% This function is workaround for that.
% Note, hard coded to round to 2 decimal places
function [X] = MyRound(x)
	X = round(x * 100) / 100;
end