% Create a cylindrical grid as in the MRST examples.
function Grid = CircleGrid(radius, radiusLayers)
    arguments
        radius = 1;
        radiusLayers = 10;
    end
    P = [];
    for r = linspace(radius/radiusLayers, radius, radiusLayers)
        [x, y, ~] = cylinder(r);
        P = [P [x(1,:); y(1,:)]];
    end
    P = unique(P', 'rows');
    Grid = pebi(triangleGrid(P));
    Grid = computeGeometry(Grid);
end