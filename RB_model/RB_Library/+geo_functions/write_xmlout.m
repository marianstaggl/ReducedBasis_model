function write_xmlout(geo_lines, file_name)
%WRITE_XMLOUT export the geo_lines as xml file

% get the selectors for hub/shroud and blades
h_sel = ismember({geo_lines.name}, 'hub');
s_sel = ismember({geo_lines.name}, 'shroud');
b_sel = find(contains({geo_lines.name}, 'blade'));

% get the geoname shortcut
[~,geo_name,~] = fileparts(file_name);

% initialize the xml file
doc_attrib.Name = geo_name; doc_attrib.BladeNumber = '0';
doc_node = xml_functions.init_xml('AiGeomInput', doc_attrib);
root_node = doc_node.getDocumentElement;

% add the nodes for the hub contour
hub_attrib.iMax = num2str(size(geo_lines(h_sel).pos,1));
xml_functions.add_node(doc_node, root_node, 'HubContour',...
    hub_attrib, geo_lines(h_sel).pos);

% add the nodes for the shroud contour
shr_attrib.iMax = num2str(size(geo_lines(s_sel).pos,1));
xml_functions.add_node(doc_node, root_node, 'ShroudContour',...
    shr_attrib, geo_lines(s_sel).pos);

% add the nodes for the blades
blade_attrib.Name = 'MainBlade'; blade_attrib.kMax = num2str(length(b_sel));
blade_attrib.Type = 'TypeCircle'; blade_attrib.ClockWise = '2';
blade_attrib.SwapHubShroud = '0';
blade_node = xml_functions.add_node(doc_node, root_node, 'Blade', ...
    blade_attrib, []);

% loop over the blade contours and add them to the xml
for i=1:length(b_sel)
    cb_attrib.iMax = num2str(size(geo_lines(b_sel(i)).pos,1));
    c_name = ['BladeProfile', num2str(i)];
    xml_functions.add_node(doc_node, blade_node, c_name,...
        cb_attrib, geo_lines(b_sel(i)).pos);
end

% write the result of the xml file
xmlwrite(file_name, doc_node)
end

