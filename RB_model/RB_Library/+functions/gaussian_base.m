function G = gaussian_base(n, lvl)
%GAUSSIAN_BASE get a base composed of shifted gaussians
% preallocate the Base G
G = zeros(n, sum(lvl)); s_vec = linspace(0,1,n); c = 1;

% set the coulumns for every level
for i=1:length(lvl)
    % get the positions of the mean
    me = linspace(0,1,lvl(i));
    sd = mean(diff(me)); if isnan(sd), sd = me; end
    
    % loop over the first level
    for j=1:lvl(i)
        % fill the column with a gaussian
        G(:,c) = gaussmf(s_vec, [sd, me(j)]); c = c + 1;
    end
end
end

