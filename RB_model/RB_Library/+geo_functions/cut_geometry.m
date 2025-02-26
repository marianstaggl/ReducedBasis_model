function [geo_lines, cutoff] = cut_geometry(geo_lines, cutoff, resamp)
%NORM_GEO cut and normalize the size of the geometry
% sort the cutoff in ascendent order
h_cut = sort(cutoff(1,:));
s_cut = sort(cutoff(2,:));

% cut hub and shroud lines
h_sel = ismember({geo_lines.name}, 'hub'); hub = geo_lines(h_sel);
s_sel = ismember({geo_lines.name}, 'shroud'); shr = geo_lines(s_sel);

% cut the hub line
[~,sel] = unique(hub.pos(:,1));
h_int = interp1(hub.pos(sel,1), hub.pos(sel,2), h_cut);
c_h = hub.pos(:,1) > h_cut(1) & hub.pos(:,1) < h_cut(2);
geo_lines(h_sel).pos = [[h_cut(1), h_int(1)];...
    hub.pos(c_h,:); [h_cut(2), h_int(2)]];

% cut the shr line
[~,sel] = unique(shr.pos(:,1));
s_int = interp1(shr.pos(sel,1), shr.pos(sel,2), s_cut);
c_s = shr.pos(:,1) > s_cut(1) & shr.pos(:,1) < s_cut(2);
geo_lines(s_sel).pos = [[s_cut(1), s_int(1)];...
    shr.pos(c_s,:); [s_cut(2), s_int(2)]];

% resample if needed
if resamp
    geo_lines(h_sel).pos = resamp_line(geo_lines(h_sel).pos, 300);
    geo_lines(s_sel).pos = resamp_line(geo_lines(s_sel).pos, 300);
end
end

function pos_new = resamp_line(pos, n)
% get a s vector for the contour line
s_old = cumsum(vecnorm(pos(2:end,:) - pos(1:end-1,:),2,2));
s_old = [0; s_old/s_old(end)]; s_new = linspace(0,1,n);

% interpolate x y and z values
pos_new(:,1) = reshape(makima(s_old, pos(:,1), s_new),[],1);
pos_new(:,2) = reshape(makima(s_old, pos(:,2), s_new),[],1);
end