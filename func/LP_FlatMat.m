%% LP_FlatMat: 
function [Vector] = LP_FlatMat(Mat)
	if size(Mat,1) ~= size(Mat,2)
		fprintf(1, 'Matrix must be symmetrical!!\n');
	else

		numRows = size(Mat,1);

		% remove lower triangle
		Vector = triu(Mat,1);

		% Generate mask
		x = ones(numRows,numRows);
		VectorRemove = tril(x);

		% reshape to vector
		Vector = reshape(Vector,numRows * numRows,1);
		VectorRemove = reshape(VectorRemove,numRows * numRows,1);
		
		% Remove lower triangle elements and diagonals from vector
		Vector(VectorRemove == 1) = [];
	end
end