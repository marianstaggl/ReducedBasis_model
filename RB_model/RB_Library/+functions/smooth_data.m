function data = smooth_data(data,std,mp)
%SMOOTH_DATA use a gaussian kernel to smooth the data

% create a s vector for the data
sv = linspace(0,1,size(data,1));
std = interp1([0, mp(1), mp(2), 1], [1e-20, std, std, 1e-20], sv);

% setup a matrix for the smoothing
G = zeros(size(data,1));

% fill the base with gaussians
for i=1:size(data,1), G(i,:) = gaussmf(sv, [std(i), sv(i)]); end

% normalize the columns
G = G./sum(G,2);

% smooth the data
for i=1:size(data,2), data(:,i) = G*data(:,i); end
end

