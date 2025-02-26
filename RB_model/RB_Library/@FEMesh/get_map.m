function [weights, dbdx_int, dbdy_int, d2bdx2_int, d2bdy2_int, inv_jac,...
    dbdx_ref, dbdy_ref] = get_map(self, varargin)
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
weights = zeros(n_elem, n_intp); inv_jac = zeros(2, 2, n_intp, n_elem);
dbdx_int = zeros(n_refp, n_intp, n_elem); dbdy_int = zeros(n_refp, n_intp, n_elem); % 1st derivatives at the integration points
dbdx_ref = zeros(n_refp, n_refp, n_elem); dbdy_ref = zeros(n_refp, n_refp, n_elem); % 1st derivatives at the reference points
d2bdx2_int = zeros(n_refp, n_intp, n_elem); d2bdy2_int = zeros(n_refp, n_intp, n_elem); % 2nd derivatives at the reference points

% get the location of the reference and int points
ref_points = self.ref_elem.ref_points;
int_points = self.ref_elem.int_points;

% get the basefunction derivatives at the integration points
dbde_int = self.ref_elem.get_dbde(int_points, 1);
dbdz_int = self.ref_elem.get_dbdz(int_points, 1);
dbde_ref = self.ref_elem.get_dbde(ref_points, 1); 
dbdz_ref = self.ref_elem.get_dbdz(ref_points, 1);

% get the basefunction derivatives at the reference points
dpdx_int = self.par_elem.get_dbde(int_points, 1);
dpdy_int = self.par_elem.get_dbdz(int_points, 1);
dpdx_ref = self.par_elem.get_dbde(ref_points, 1); 
dpdy_ref = self.par_elem.get_dbdz(ref_points, 1);

% loop over all of the elements
for i=1:length(elem_id)
    c_elem = elem_id(i);
    
    % get the coordinates of the support nodes
    nodes = self.nodes(self.elems(c_elem,:),:);
    
    % map the derivatives onto the deformed element
    [jac_det, dbdx_int(:,:,c_elem), dbdy_int(:,:,c_elem),...
        d2bdx2_int(:,:,c_elem), d2bdy2_int(:,:,c_elem),...
        inv_jac(:,:,:,c_elem), dbdx_ref(:,:,c_elem), dbdy_ref(:,:,c_elem)] = ...
        map_single_elem(n_intp, n_refp, self.ref_elem.order, dpdx_int, dpdy_int,...
        dpdx_ref, dpdy_ref, dbde_int, dbdz_int, dbde_ref, dbdz_ref, nodes );
    
    % multiply the jacobian determinant by its nodal weights
    weights(c_elem,:) = jac_det.*self.ref_elem.int_weights;
end
end

function  [jac_det, dbdx_int, dbdy_int, d2bdx2_int, d2bdy2_int, inv_jac_int,...
    dbdx_ref, dbdy_ref] = map_single_elem( n_intp, n_refp, order, ...
    dpdx_int, dpdy_int, dpdx_ref, dpdy_ref, dbde_int, dbdz_int, ...
    dbde_ref, dbdz_ref, support_nodes )
% map the properties of a single element

% preallocate the arrays
jac_det = zeros(1,n_intp); inv_jac_int = zeros(2, 2, n_intp);
dbdx_int = zeros(n_refp, n_intp); dbdy_int = zeros(n_refp, n_intp);
dbdx_ref = zeros(n_refp, n_refp); dbdy_ref = zeros(n_refp, n_refp);

% loop over all of the integration points
for j=1:n_intp
    % get the jacobian matrix in integration and reference points
    jac_int = (support_nodes'*[dpdx_int(:,j), dpdy_int(:,j)])';
    jac_det(j) = det(jac_int);
    
    % check if the jacobian matrix is positive
    if jac_det(j)<0, error('negative jacobian determinant! sort nodes...'); end
    
    % invert the jacobian matricies on reference and integration points
    ijac_int = inv(jac_int);
    inv_jac_int(:,:,j) = ijac_int; % save the jacobian
    
    % map the derivatives of the integration points
    dbdx_int(:,j) = ijac_int(1,1)*dbde_int(:,j) + ijac_int(1,2)*dbdz_int(:,j);
    dbdy_int(:,j) = ijac_int(2,1)*dbde_int(:,j) + ijac_int(2,2)*dbdz_int(:,j);
end

% loop over all of the integration points
for j=1:n_refp
    % get the jacobian matrix at the reference points
    jac_ref = (support_nodes'*[dpdx_ref(:,j), dpdy_ref(:,j)])';
    
    % invert the jacobian matricies on reference and integration points
    ijac_ref = inv(jac_ref); 
    
    % map the derivatives of the reference points
    dbdx_ref(:,j) = ijac_ref(1,1)*dbde_ref(:,j) + ijac_ref(1,2)*dbdz_ref(:,j);
    dbdy_ref(:,j) = ijac_ref(2,1)*dbde_ref(:,j) + ijac_ref(2,2)*dbdz_ref(:,j);
end

% return the first and second derivatives at the int points
if order < 2
    d2bdx2_int = zeros(size(dbdx_int));
    d2bdy2_int = zeros(size(dbdy_int));
else
    d2bdx2_int = (dbdx_int'*dbdx_ref')'; 
    d2bdy2_int = (dbdy_int'*dbdy_ref')'; 
end
end