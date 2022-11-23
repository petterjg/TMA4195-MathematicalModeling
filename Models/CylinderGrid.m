% Create a cylindrical grid as in the MRST examples.
function Grid = CylinderGrid(radius, height, radiusLayers, heightLayers)
    arguments
        radius = 1;
        height = 1;
        radiusLayers = 10;
        heightLayers = 10;
    end
    P = [];
    for r = linspace(radius/radiusLayers, radius, radiusLayers)
        [x, y, ~] = cylinder(r);
        P = [P [x(1,:); y(1,:)]];
    end
    P = unique(P', 'rows');
    aG = pebi(triangleGrid(P));
    Grid = makeLayeredGrid(aG, ones(heightLayers,1)*height/heightLayers);
    Grid = computeGeometry(Grid);
    Grid = createAugmentedGrid(Grid);
end