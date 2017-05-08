function [] = PlotMotion(mov,fd)

%
% This function will plot the translation and rotation parameters following
% realignment and fd
% AKA a simplified version of ThePlot without the time series data.
%
% ------
% INPUTS
% ------
% 
% mov   - an N x 6 matrix, where N = number of time points. Each column
%       represents a different motion parameter as per the following order:
%       1 - x translation; 2 - y translation; 3 - z translation; 4 - pitch;
%       5 - roll; 6 - yaw.
%
% fd	- a framewise displacement measure (e.g., from Power or Jenkinson or Van Dijk)
% -------
% OUTPUTS
% -------
%
%==========================================================================

N = size(mov,1);

figure
subplot(3,1,1)
plot(mov(:,1));
hold on
plot(mov(:,2),'g');
plot(mov(:,3),'r');
title('translation','fontsize',15,'fontweight','bold')
legend({'x','y','z'})
ylabel('mm')
xlabel('time (volumes)')
xlim([1 N])

subplot(3,1,2)
plot(mov(:,4));
hold on
plot(mov(:,5),'g');
plot(mov(:,6),'r');
title('rotation','fontsize',15,'fontweight','bold');
legend({'pitch','roll','yaw'})
ylabel('degrees')
xlabel('time (volumes)')
xlim([1 N])

subplot(3,1,3)
plot(fd);
title('fd','fontsize',15,'fontweight','bold');
ylabel('mm')
xlabel('time (volumes)')
xlim([1 N])
