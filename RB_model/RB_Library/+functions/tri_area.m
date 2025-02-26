function A = tri_area(T, x, y, z)
%TRI_AREA get the areas of triangles in a triangulated surface
%   Detailed explanation goes here

% get the nodes of the triangles
X = x(T); Y = y(T); Z = z(T);

% get the points of the triangles
P1 = [X(:,1), Y(:,1), Z(:,1)]; 
P2 = [X(:,2), Y(:,2), Z(:,2)];
P3 = [X(:,3), Y(:,3), Z(:,3)];

% get the area of all the triangles
A = 1/2*vecnorm(cross((P2-P1)',(P3-P1)'));
end

