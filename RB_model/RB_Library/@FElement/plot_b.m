function plot_b(self)
%PLOT_B plot the elements base functions

% define the number of gridpoints
num_g = 100;

% create a plotting domain and reshape it
[ETA, ZETA] = meshgrid(linspace(-1,1,num_g), linspace(-1,1,num_g));
eta = reshape(ETA,[],1); zeta = reshape(ZETA,[],1);

% get the correspoinding functional values
b = self.get_b([eta,zeta]);

% loop over the basefunctions
for i=1:size(b,1)
    
    % reshape the basefunctions
    B = reshape(b(i,:), num_g, num_g);
    
    % plot the corresponding basefunction
    subplot(self.order+1, self.order+1, i), contourf(ETA,ZETA,B,20); colorbar()
    title([num2str(i) '. basefunction']); xlabel('eta'); ylabel('zeta'); axis equal
end
end

