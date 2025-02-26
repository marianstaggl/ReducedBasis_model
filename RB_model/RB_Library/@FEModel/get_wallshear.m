function wall_shear = get_wallshear(self, curr_cond, opts)
%GET_Y_PLUS calculate y plus for each meshpoints

% get the vectors with x and y velocity
curr_cond = reshape(curr_cond, [], numel(self.sys_variables));
uid = ismember(self.sys_variables,'u'); uvec = curr_cond(:,uid);
vid = ismember(self.sys_variables,'v'); vvec = curr_cond(:,vid);
wall_shear = zeros(size(curr_cond,1),1);
id_count = zeros(size(curr_cond,1),1);

% loop over all the lines
for i=1:numel(self.mesh.lines)
    % select one line from the mesh
    line = self.mesh.lines(i);
    n1 = line.n1; n2 = line.n2;
    elems = line.elems;
    subs = line.sub_id;
    
    % loop over all the elements
    for j=1:length(elems)
        % get the current element and its mappings
        curr_elem = self.mesh.elems(elems(j),:);
        dbdx = self.dbdx(:,subs(:,j),elems(j));
        dbdy = self.dbdy(:,subs(:,j),elems(j));
        uval = uvec(curr_elem); vval = vvec(curr_elem);
        
        % get the current velocitys at the quad points
        dudx = sum(dbdx.*uval,1); dudy = sum(dbdy.*uval,1);
        dvdx = sum(dbdx.*vval,1); dvdy = sum(dbdy.*vval,1);
        
        % get the wall normal velocity derivative
        shear = opts.mue*abs((dudx*n2(j,1) + dvdx*n2(j,2))*n1(j,1) + ...
            (dudy*n2(j,1) + dvdy*n2(j,2))*n1(j,2));
        
        % save the wall shear
        c_ids = curr_elem(subs(:,j));
        wall_shear(c_ids) = wall_shear(c_ids) + shear';
        id_count(c_ids) = id_count(c_ids) + 1;
    end
end

% average double values
sel = id_count > 0;
wall_shear(sel) = wall_shear(sel) ./ id_count(sel);
end

