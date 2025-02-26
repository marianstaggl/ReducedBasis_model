function [sel, wei] = slice_tri(T, X, sv, sp)
% SLICE_TRI slice a triangulated surface along x or y axis
%   Input values:
%   T: triangulation of the surface
%   X: coordinates of the points in a nx3 array
%   sv: translation of the slicing plane
%   sp: 1 or 2 (defines the direction of the slice)
%
%   Output values:
%   sel: the nodes making up the edges with an intersection
%   wei: weighting values for the nodes in sel for interpolation

% cut the points into 2 parts
y_raw = X(:,sp) - sv; y_tri = y_raw(T);

% find the intersection of all possible edges
e12 = y_tri(:,1)./(y_tri(:,1) - y_tri(:,2)); int12 = (e12 >= 0) & (e12 <= 1);
e23 = y_tri(:,2)./(y_tri(:,2) - y_tri(:,3)); int23 = (e23 >= 0) & (e23 <= 1);
e31 = y_tri(:,3)./(y_tri(:,3) - y_tri(:,1)); int31 = (e31 >= 0) & (e31 <= 1);

% find the edges which intersect the r_threshold
sel1 = T(int12,[1,2]); weig1 = [1-e12(int12), e12(int12)];
sel2 = T(int23,[2,3]); weig2 = [1-e23(int23), e23(int23)];
sel3 = T(int31,[3,1]); weig3 = [1-e31(int31), e31(int31)];

% get them in one array (caution, there are duplicate edges)
sel = [sel1; sel2; sel3]; wei = [weig1; weig2; weig3];

% check if there are nans in the weights
if any(isnan(wei)), error('nan in weightvec, check again...'); end
end