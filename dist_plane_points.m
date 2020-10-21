function [dist, signed_dist] = dist_plane_points(...
    plane_origin, plane_vect1, plane_vect2, points)
% [dist, signed_dist] = dist_plane_points(...
%   plane_origin, plane_vect1, plane_vect2, points)
% compute distances between plane and points
%   plane_origin  - 1*3 array - one point      that belongs to the plane
%   plane_vect1   - 1*3 array - one vector     that belongs to the plane
%   plane_vect2   - 1*3 array - another vector that belongs to the plane
%   points        - n*3 array - points to compute distances to
%   dist          - n*1 array - copmuted distances
%   signed_dist   - n*1 array - copmuted distances with sign: 
%                               positive away from zero of coordinate system,
%                               negative towards   zero of coordinate system

% compute a vector product - a normal vector to the plane
abc = cross(plane_vect1, plane_vect2);
if norm(abc) < 1e-10 .* norm([ plane_vect1, plane_vect2 ])
    warning('dist_plane_points: plane_vect1 and plane_vect2 are almost collinear, plane is defined poorly')
end

% find plane equation in form of a*x+b*y+c*z = p
p = sum(plane_origin .* abc);

% normalize the plane equation 
k = norm(abc);
if (p<0)
    k = -k;
end
abc = abc / k;
p = p / k;

% check
%d1 = sum(abc .* plane_origin) - p
%d2 = sum(abc .* (plane_origin + plane_vect1)) - p
%d3 = sum(abc .* (plane_origin + plane_vect2)) - p

% calc distance
signed_dist = sum(repmat(abc, size(points, 1), 1) .* points, 2) - p;
dist = abs(signed_dist);



