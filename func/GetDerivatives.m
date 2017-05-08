%% GetDerivatives: function description
function [out] = GetDerivatives(ts,getSquares)

     % This function computes the temporal derivates (i.e., backwards differences)
     % as well as the squares of a set of nuissance regressors
     %
     % Linden Parkes, Brain & Mental Health Laboratory, 2016
     % ------------------------------------------------------------------------------

     if nargin < 2
          getSquares = 1;
     end

     % detrend
     ts_detr = detrend(ts,'linear');

     % temporal derivatives
     ts_diff = [
               zeros(1,size(ts_detr,2));
               diff(ts_detr)
               ];

     % concatenate
     x = [ts_detr ts_diff];

     if getSquares == 1
          % squared
          x_sq = x.^2;

          % concatenate
          out = [x x_sq];
     else
          out = x;
     end
          
end