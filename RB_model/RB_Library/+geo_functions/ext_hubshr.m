function geo_lines = ext_hubshr(geo_lines, ext_f, order, np)
%EXT_HUBSHR extend hub and shroud line at the backend

% select hub and shroudline from the geo
hsel = ismember({geo_lines.name}, 'hub'); hub = geo_lines(hsel);
ssel = ismember({geo_lines.name}, 'shroud'); shr = geo_lines(ssel);

% transfer the ext into x coordinates
new_lim = ext_f*(hub.pos(end,1) - hub.pos(1,1)) + hub.pos(end,1);

% get the last np points used for extrapolation
ph = polyfit(hub.pos(end-np:end,1), hub.pos(end-np:end,2), order);
ps = polyfit(shr.pos(end-np:end,1), shr.pos(end-np:end,2), order);

% extend the x of hub and shroud
xh = hub.pos(end,1):median(diff(hub.pos(:,1))):new_lim;
xs = shr.pos(end,1):median(diff(shr.pos(:,1))):new_lim;
yh = polyval(ph, xh); yh = (yh(2:end) - yh(1)) + hub.pos(end,2);
ys = polyval(ps, xs); ys = (ys(2:end) - ys(1)) + shr.pos(end,2);

% extend the hub and shroudline
geo_lines(hsel).pos = [geo_lines(hsel).pos; [xh(2:end)', yh']];
geo_lines(ssel).pos = [geo_lines(ssel).pos; [xs(2:end)', ys']];
end

