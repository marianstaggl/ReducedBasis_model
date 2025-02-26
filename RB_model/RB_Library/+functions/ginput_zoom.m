function [x,y] = ginput_zoom()
%GINPUT_ZOOM ginput function with zoom capability

% preallocate the x and y value
x = []; y = [];

% display the text for the point selection
text(0.1, 0.7, {'select points', 'reset zoom: enter',...
    'break: del + enter'}, 'units', 'normalized');

% start a while loop
while true
    % allow zoom to mark the points
    zoom on; pause(); zoom off;
    
    % use ginput to select points
    [x_tmp,y_tmp,b_tmp] = ginput();
    
    % save the newly added points
    x = [x; x_tmp(b_tmp~=8)]; 
    y = [y; y_tmp(b_tmp~=8)];
    
    % check if b is empty
    if isempty(b_tmp), continue; end
    
    % check the key pressed
    if b_tmp(end) == 8, break; end
end

end

