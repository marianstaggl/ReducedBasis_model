function [M, Bl, Bu, G, ids, names] = inletbc_base(n)
%INLETBC_BASE Create a base of profiles for the inlet conditions

% get the mean and lean profiles
M = mean_prof(n);

% get the boundarylayer profiles
B = bound_prof(n); Bl = B; Bu = flipud(B);

% set the G base to be a base of gaussian functions
G = gauss_prof(n);

% define an id matrix for the different blocks
ids = [0 cumsum([size(M,2), size(Bl,2), size(Bu,2), size(G,2)])];

% set the names of the columns
names = setup_columnnames(M, Bl, Bu, G);
end

function names = setup_columnnames(~, Bl, Bu, G)
% get the names of the M parameters
m_names = {'mean_flow', 'lean_flow'};

% get the names of the bl parameters
for i=1:size(Bl,2), bl_names{i} = ['hub_blayer_' num2str(i/size(Bl,2))]; end

% get the names of the bu parameters
for i=1:size(Bu,2), bu_names{i} = ['shr_blayer_' num2str(i/size(Bu,2))]; end

% get the name of the localized parameters
for i=1:size(G,2), g_names{i} = ['localized_' num2str(i/size(G,2))]; end

% merge them together
names = horzcat(m_names, bl_names, bu_names, g_names);

% replace the points with underscores for excel output
names = strrep(names, '.', '_');
end

function G = gauss_prof(n)

% approximate the boundary region with gaussians
s_min = 0.025; s_max = 1-0.025; s_vec = linspace(0, 1, n); 
m_vec = linspace(s_min, s_max, 32); sdl = mean(diff(m_vec));

% set the coulumns for every level
for i=1:length(m_vec), G(:,i) = gaussmf(s_vec, [sdl, m_vec(i)]); end
end

function M = mean_prof(n)
% get a profile to encode the mean/lean
M = [ones(n,1), linspace(-1,1,n)'];

% normalize the profiles by max values
M = M./max(M,1);
end

% function B = bound_prof(n)
% 
% % approximate the boundary region with gaussians
% s_min = 0; s_max = 0.1; s_vec = linspace(0, 1, n); 
% m_vec = linspace(s_min, s_max, 10); sdl = mean(diff(m_vec));
% 
% % set the coulumns for every level
% for i=1:length(m_vec), B(:,i) = gaussmf(s_vec, [sdl, m_vec(i)]); end
% end

function B = bound_prof(n)
% approximate the boundary region with gaussians
sd_min = 1e-3; sd_max = 0.035; s_vec = linspace(0,1,n);
sdl = (sd_max - sd_min)*(1-cos(linspace(0, pi/2, 8))) + sd_min;

% set the coulumns for every level
for i=1:length(sdl), B(:,i) = gaussmf(s_vec, [sdl(i), 0]); end
end

% function B = bound_prof(n)
% % get a power base to encode the boundary thickness
% [~,B] = bc_functions.power_law(n, [8, 13, 21, 34]); B = 1-B;
% end

% function B = bound_prof(n)
% % get a nondimensionalized boundary layer (wall law)
% [y_plus,u_plus] = bc_functions.wall_law(1e3);
% 
% % normalize y plus and u plus by max values
% y_plus = y_plus/max(y_plus); u_plus = u_plus/max(u_plus);
% 
% % define the b layer thickness levels
% n_lvl = 5; min_th = 1e-5; max_th = 0.2; s_vec = linspace(0,1,n);
% th_lvl = (max_th - min_th)*(1-cos(linspace(0, pi/2, n_lvl))) + min_th;
% 
% % loop over the levels
% for i=1:length(th_lvl)
%     % adjust the boundary thickness
%     y_tmp = [th_lvl(i)*y_plus, 1]; u_tmp = [u_plus, 1];
%     
%     % interpolate at the span vector
%     B(:,i) = 1 - interp1(y_tmp, u_tmp, s_vec);
% end
% end
