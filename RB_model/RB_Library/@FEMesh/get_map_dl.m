function [weights, dbdx_int, dbdy_int, d2bdx2_int, d2bdy2_int, inv_jac,...
    dbdx_ref, dbdy_ref] = get_map_dl(self, varargin)
%GET_MAP map the derivatives and the integration weights
%   to evaluate the inner products of the residual we need to perform a
%   numerical integration. as the integrational weights as well as the
%   derivatives are only defined at the reference element we need to map
%   them onto the new element. we map the first as well as the second
%   derivative.

% if no elements are specified, loop over all of them
if isempty(varargin), elem_id = 1:1:self.n_elem;
else, elem_id = varargin{1}; end

% create shortcuts to the number of int and ref points
n_elem = self.n_elem;
n_intp = size(self.ref_elem.int_points,1);
n_refp = size(self.ref_elem.ref_points,1);

% preallocate the arrays
weights = zeros(n_elem, n_intp);
inv_jac = zeros(2, 2, n_intp, n_elem);
dbdx_int = zeros(n_refp, n_intp, n_elem); 
dbdy_int = zeros(n_refp, n_intp, n_elem); % 1st derivatives at the integration points
dbdx_ref = zeros(n_refp, n_refp, n_elem); 
dbdy_ref = zeros(n_refp, n_refp, n_elem); % 1st derivatives at the reference points
d2bdx2_int = zeros(n_refp, n_intp, n_elem); 
d2bdy2_int = zeros(n_refp, n_intp, n_elem); % 2nd derivatives at the reference points

% check if the nodes are dlarrays, if so, convert the rest
if isa(self.nodes, 'dlarray')
    weights = dlarray(weights); inv_jac = dlarray(inv_jac); 
    dbdx_int = dlarray(dbdx_int); dbdy_int = dlarray(dbdy_int);
    dbdx_ref = dlarray(dbdx_ref); dbdy_ref = dlarray(dbdy_ref);
    d2bdx2_int = dlarray(d2bdx2_int); d2bdy2_int = dlarray(d2bdy2_int);
end

% get the basefunction derivatives at the integration points
dbde_int = self.ref_elem.get_dbde(self.ref_elem.int_points, 1);
dbdz_int = self.ref_elem.get_dbdz(self.ref_elem.int_points, 1);
dbde_ref = self.ref_elem.get_dbde(self.ref_elem.ref_points, 1); 
dbdz_ref = self.ref_elem.get_dbdz(self.ref_elem.ref_points, 1);

% get the basefunction derivatives at the reference points
dpdx_int = self.par_elem.get_dbde(self.ref_elem.int_points, 1);
dpdy_int = self.par_elem.get_dbdz(self.ref_elem.int_points, 1);
dpdx_ref = self.par_elem.get_dbde(self.ref_elem.ref_points, 1); 
dpdy_ref = self.par_elem.get_dbdz(self.ref_elem.ref_points, 1);

% loop over all of the elements
for i=1:length(elem_id)
    c_elem = elem_id(i);
    
    % get the coordinates of the support nodes
    nodes = self.nodes(self.elems(c_elem,:),:);
    
    % map the derivatives onto the deformed element
    [jac_det, dbdx_int(:,:,c_elem), dbdy_int(:,:,c_elem),...
        d2bdx2_int(:,:,c_elem), d2bdy2_int(:,:,c_elem),...
        inv_jac(:,:,:,c_elem), dbdx_ref(:,:,c_elem), dbdy_ref(:,:,c_elem)] = ...
        map_single_elem( self.ref_elem.order, dpdx_int, dpdy_int,...
        dpdx_ref, dpdy_ref, dbde_int, dbdz_int, dbde_ref, dbdz_ref, nodes );
    
    % multiply the jacobian determinant by its nodal weights
    weights(c_elem,:) = jac_det.*self.ref_elem.int_weights;
end
end

function  [jac_det, dbdx_int, dbdy_int, d2bdx2_int, d2bdy2_int, inv_jac_int,...
    dbdx_ref, dbdy_ref] = map_single_elem( order, dpdx_int, dpdy_int, dpdx_ref,...
    dpdy_ref, dbde_int, dbdz_int, dbde_ref, dbdz_ref, support_nodes )

% reshape the arrays such that the integration point is in 3rd dimension
dp_int = cat(2, [permute(dpdx_int, [1, 3, 2]),...
    permute(dpdy_int, [1, 3, 2])]);
dp_ref = cat(2, [permute(dpdx_ref, [1, 3, 2]),...
    permute(dpdy_ref, [1, 3, 2])]);

% calculate the jacobian in every integration point
jac_int = pagemtimes(repmat(support_nodes', 1, 1, 4), dp_int);
jac_int_det = jac_int(1,1,:).*jac_int(2,2,:) - ...
    jac_int(1,2,:).*jac_int(2,1,:);
inv_jac_int = 1./jac_int_det .* [jac_int(2,2,:), - jac_int(2,1,:); ...
    - jac_int(1,2,:) , jac_int(1,1,:)];
dbdx_int = reshape(inv_jac_int(1,1,:), 1, []) .* dbde_int + ...
    reshape(inv_jac_int(1,2,:), 1, []) .* dbdz_int;
dbdy_int = reshape(inv_jac_int(2,1,:), 1, []) .* dbde_int + ...
    reshape(inv_jac_int(2,2,:), 1, []) .* dbdz_int;

% calculate the mapping in every reference point
jac_ref = pagemtimes(repmat(support_nodes', 1, 1, 4), dp_ref);
jac_ref_det = jac_ref(1,1,:).*jac_ref(2,2,:) - ...
    jac_ref(1,2,:).*jac_ref(2,1,:);
inv_jac_ref = 1./jac_ref_det .* [jac_ref(2,2,:), - jac_ref(2,1,:); ...
    - jac_ref(1,2,:) , jac_ref(1,1,:)];
dbdx_ref = reshape(inv_jac_ref(1,1,:), 1, []) .* dbde_ref + ...
    reshape(inv_jac_ref(1,2,:), 1, []) .* dbdz_ref;
dbdy_ref = reshape(inv_jac_ref(2,1,:), 1, []) .* dbde_ref + ...
    reshape(inv_jac_ref(2,2,:), 1, []) .* dbdz_ref;

% calculate the values at the reference points
if any(jac_int_det < 0) || any(jac_ref_det < 0)
    error('negative jacobian determinant! sort nodes...'); 
end

% reshape the jacobian determinant
jac_det = reshape(jac_int_det, 1, []);

% return the first and second derivatives at the int points
if order < 2
    d2bdx2_int = zeros(size(dbdx_int));
    d2bdy2_int = zeros(size(dbdy_int));
else
    d2bdx2_int = (dbdx_int'*dbdx_ref')'; 
    d2bdy2_int = (dbdy_int'*dbdy_ref')'; 
end
end