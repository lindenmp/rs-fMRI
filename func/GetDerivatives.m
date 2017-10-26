%% GetDerivatives: function description
function [out] = GetDerivatives(ts,detr,getSquares)

     % This function computes the temporal derivates (i.e., backwards differences)
     % as well as the squares of a set of nuissance regressors
     %
     % Linden Parkes, Brain & Mental Health Laboratory, 2016
     % ------------------------------------------------------------------------------
     
     if nargin < 2
          detr = 1;
     end

     if nargin < 3
          getSquares = 1;
     end

     % detrend
     if detr == 1
          ts = detrend(ts);
     end

     % temporal derivatives
     ts_diff = [
               zeros(1,size(ts,2));
               diff(ts)
               ];

     % concatenate
     x = [ts ts_diff];

     if getSquares == 1
          % squared
          x_sq = x.^2;

          % concatenate
          out = [x x_sq];
     else
          out = x;
     end
          
end