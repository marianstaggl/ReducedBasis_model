function [vals, names] = par2profile(base, par)
%PAR2PROFILE this function sets up an artificial velocity profile at the
%inlet of the TCF and returns the axial, circumferential and radial vels
%   input:
%   par --> parameter struct containing the fields: 'vel_mean', 'vel_lean',...
%   n --> number of lines for the profiles

% preallocate the vals cell
vals = cell(1,numel(par));
names = cell(1,numel(par));

% loop over the parameters
for i=1:numel(par)
    vals{i} = base*(par(i).data)';
    names{i} = par(i).name;
end

% search for vel and angles and convert them to cax, crad and ccirc
[vals, names] = convert2velcomp(vals, names);
end

function [vals, names] = convert2velcomp(vals, names)

% check if there are the columns vel, yaw and pitch
[~,ids] = ismember({'vel', 'yaw', 'pitch'}, names);

%  if they are not available return
if any(~ids), return; end

% get the components vel yaw and pitch
vel = vals{ids(1)}; yaw = vals{ids(2)}; pitch = vals{ids(3)};

% otherwise transform to vel components
cf = 1./sqrt(1 + tand(yaw).^2 + tand(pitch).^2);
vals{end+1} = vel.*cf; vals{end+1} = vel.*tand(yaw); vals{end+1} = vel.*tand(pitch);
names = horzcat(names, {'c_ax', 'c_circ', 'c_rad'});
end

% 
% % preallocate the velocity-component arrays
% vel = zeros(n, numel(par), 3); 
% yaw = zeros(n, numel(par), 3); 
% pitch = zeros(n, numel(par), 3);
% 
% % loop over the parametrization and get the first
% for i=1:numel(par)
%     % get the first accuracy level of the profiles
%     [vel(:,i,1), yaw(:,i,1), pitch(:,i,1)] = get_1st_level(par(i), n);
%     
%     % get the second accuracy level of the profiles
%     [vel(:,i,2), yaw(:,i,2), pitch(:,i,2)] = get_2nd_level(par(i), n);
%     
%     % get the third accuracy level of the profiles
%     [vel(:,i,3), yaw(:,i,3), pitch(:,i,3)] = get_3rd_level(par(i), n);
% end
% 
% % sum up the different flow levels up until the desired level
% vel = sum(vel(:,:,1:ip.Results.acc_levels(1)),3);
% yaw = sum(yaw(:,:,1:ip.Results.acc_levels(2)),3);
% pitch = sum(pitch(:,:,1:ip.Results.acc_levels(3)),3);
% 
% % transform the profiles into velocity components
% cf = 1./sqrt(1 + tand(yaw).^2 + tand(pitch).^2);
% c_ax = vel.*cf; c_circ = vel.*tand(pitch); c_rad = vel.*tand(yaw);
% 
% % names of the columns
% names = {'c_ax', 'c_circ', 'c_rad'};
% vals = {c_ax, c_circ, c_rad};
% end
% 
% function [vel, yaw, pitch] = get_1st_level(par, n)
% % get the bases for the mean flow features
% M = bc_functions.inletbc_base(n);
% 
% % get the profiles for yaw, pitch and velocity
% yaw = M*[par.yaw_mean; par.yaw_lean];
% pitch = M*[par.pitch_mean; par.pitch_lean];
% vel = M*[par.vel_mean; par.vel_lean];
% end
% 
% function [vel, yaw, pitch] = get_2nd_level(par, n)
% % get the bases for the boundary flow features
% [~,Bl,Bu] = bc_functions.inletbc_base(n);
% B_base = [Bl, Bu];
% 
% % get the profiles for yaw, pitch and velocity
% yaw = B_base*[reshape(par.lo_yawbl,[],1); reshape(par.up_yawbl,[],1)];
% pitch = B_base*[reshape(par.lo_pitchbl,[],1); reshape(par.up_pitchbl,[],1)];
% vel = B_base*[reshape(par.lo_velbl,[],1); reshape(par.up_velbl,[],1)];
% end
% 
% function [vel, yaw, pitch] = get_3rd_level(par, n)
% % get teh bases for the localized flow features
% [~,~,~,G] = bc_functions.inletbc_base(n);
% 
% % get the profiles for yaw, pitch and velocity
% yaw = G*reshape(par.yaw_lf,[],1);
% pitch = G*reshape(par.pitch_lf,[],1);
% vel = G*reshape(par.vel_lf,[],1);
% end