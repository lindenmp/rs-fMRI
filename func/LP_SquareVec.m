%% LP_SquareVec: 
function [Mat] = LP_SquareVec(Vector,numRowCol)

	Mat = triu(ones(numRowCol),1);
	Mat(Mat == 1) = Vector;

end