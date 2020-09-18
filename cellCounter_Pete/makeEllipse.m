function [xNew, yNew] = makeEllipse(a, b, center, theta, vertPoints)
    % subroutine to update x & y vertices values of the ellipse from
    % trig parametric equation (see Wikipedia for equation details)
    xNew = center(1) + a*cosd(theta)*cosd(vertPoints) - b*sind(theta)*sind(vertPoints); % x-coords of vertices
    yNew = center(2) + a*sind(theta)*cosd(vertPoints) + b*cosd(theta)*sind(vertPoints); % y-coords
end