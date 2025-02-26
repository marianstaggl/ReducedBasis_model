function [outputArg1,outputArg2] = get_turb_BaldwinLomax(self, curr_cond, opts)
%GET_TURB_BALDWINLOMAX use a baldwin lomax model to calculate the turb visc
%   The model uses y+ and the vorticity to calculate the turbulent
%   viscosity

% define the constants for the model (Source:
% https://www.cfd-online.com/Wiki/Baldwin-Lomax_model)
A_plus = 26; C_cp = 1.6; C_kleb = 0.3;
C_wk = 0.25; k = 0.4; K = 0.0168;
rho_ref = 1.1841; p_ref = 101325;

% preallocate some variables
curr_cond = reshape(curr_cond, [], 3); n_elem = size(self.mesh.elems, 1);
uid = ismember(self.sys_variables,'u'); uvec = curr_cond(:,uid);
vid = ismember(self.sys_variables,'v'); vvec = curr_cond(:,vid);
aid = ismember(self.sys_variables,'a'); avec = curr_cond(:,aid);

% get the basefunction at the integration points
int_points = self.mesh.ref_elem.int_points;
b = self.mesh.ref_elem.get_b(int_points);
vel = sqrt(uvec.^2 + vvec.^2); u_diff = max(vel) - min(vel);

% calculate y+ for each node in the domain
yvec = self.wall_dist; y_max = max(yvec);
wall_shear = self.get_wallshear(curr_cond, opts);
ypvec = yvec .* sqrt(wall_shear) / opts.mue;

% loop over the elements
for i=1:n_elem
    % get the current element and its mappings
    curr_elem = self.mesh.elems(i,:);
    
    % get the current velocitys at the quad points
    dudy = sum(self.dbdy(:,:,i).*uvec(curr_elem),1);
    dvdx = sum(self.dbdx(:,:,i).*vvec(curr_elem),1);
    y = sum(b.*yvec(curr_elem),1); yp = sum(b.*ypvec(curr_elem),1);
    a = sum(b.*avec(curr_elem),1); 
    
    % calculate the density
    rho = (a.^2/(gamma*p_ref)*(rho_ref.^(gamma))).^(1/(gamma-1));
    
    %% calculate mue_t_inner
    % get the vorticity at the quad points
    omega_12 = 1/2*(dudy - dvdx); omega_21 = 1/2*(dvdx - dudy);
    omega_n = sqrt(2*(omega_12.^2 + omega_21.^2));
    
    % get the y+ value at the quad points
    len = k*y.*(1 - exp(-yp/A_plus));
    
    % calculate the inner mue turb
    i_mue_t = rho.*len .^2.*omega_n;
    
    %% calculate mue_t_outer
    % get the f wake function
    F_max = y_max * omega_n * (1 - exp(-yp(A_plus)));
    F_wake = min([y_max*F_max, C_wk*y_max*u_diff^2/F_max]);
    F_kleb = (1 + 5.5*(y*C_kleb/y_max)^6)^(-1);
    
    % calculate the outer mue turb
    o_mue_t = K*C_cp*rho.*F_wake.*F_kleb;
end

% save the turbulence field in the options
opts.mue_t = mue_t
end

