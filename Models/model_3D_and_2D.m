close all

mrstModule add ad-core mrst-gui 

jsonfile = fileread('input.json');
jsonstruct = jsondecode(jsonfile);

paramobj = ReactionDiffusionInputParams_new(jsonstruct);


% Define grid
radius = 0.22*micro*meter;
height = 15*nano*meter;
radiusLayers = 10;
heightLayers = 10;

% Initialize variables
avo = 6.022e23;
dim = 3; % Choose between the 3D and 2D model here!
if dim == 3
    % 3D
    Grid = CylinderGrid(radius,height,radiusLayers,heightLayers);
    vols = Grid.cells.volumes;
    nc = Grid.cells.num;
    topLayerIndices = (1 : nc/heightLayers);
    bottomLayerIndices = (nc - nc/heightLayers + 1 : nc);
    cR0 = 1000/(avo*micro^2)/(height/heightLayers);
end
if dim == 2
    % 2D
    Grid = CircleGrid(radius,radiusLayers);
    vols = Grid.cells.volumes*height;
    nc = Grid.cells.num;
    topLayerIndices = (1 : nc);
    bottomLayerIndices = topLayerIndices;
    cR0 = 1000/(avo*micro^2)/(height);
end

% Spread the 5000 ns across the cells in the center of the top layer.
cN0 = 5000/(avo*sum(vols(find(vols(topLayerIndices) < min(vols)*1.01))));



paramobj.k_on = 4e3;
paramobj.k_off = 5;
paramobj.N.D = 3e-10;
paramobj.R.D = 0;
paramobj.RN.D = 0;
paramobj.G = Grid;

paramobj = paramobj.validateInputParams_new();

model = ReactionDiffusion_new(paramobj);


% Setup schedule
total = 1e-4;
n  = 2000; % Reduced for faster runtime.
dt = total/n;
step = struct('val', dt*ones(n, 1), 'control', ones(n, 1));

control.none = [];
schedule = struct('control', control, 'step', step);

% Setup initial state

initcase = 1;
switch initcase
  case 1
    cR = zeros(nc, 1);
    cR(bottomLayerIndices) = cR0;
    cN = zeros(nc, 1);
    cN(find(vols(topLayerIndices) < min(vols)*1.01)) = cN0;
    cRN = zeros(nc, 1);
  case 2
    cR = ones(nc, 1);
    cN = ones(nc, 1);
    cRN = zeros(nc, 1);
end

initstate.R.c = cR;
initstate.N.c = cN;
initstate.RN.c = cRN;

% run simulation

nls = NonLinearSolver();
nls.errorOnFailure = false;

[~, states, report] = simulateScheduleAD(initstate, model, schedule, 'NonLinearSolver', nls);

%%
% Plot 3D
% Remove empty states (could have been created if solver did not converge)
ind = cellfun(@(state) ~isempty(state), states);
states = states(ind);

figure(1); figure(2); figure(3);

for istate = 2000%1 : numel(states)

    state = states{istate};

    set(0, 'currentfigure', 1);
    cla
    view(30,-50);
    cols = state.R.c;
    plotCellData(model.G, cols);
    colorbar
    caxis([0, max(cols)]);
    title('R concentration')
    set(gca,'FontSize',14)
    exportgraphics(gcf,'Project/Figures/cR_it2000_3D.png','ContentType','vector')
    
    
    set(0, 'currentfigure', 2);
    cla
    view(30,-50);
    cols = state.N.c;
    plotCellData(model.G, cols);
    colorbar
    caxis([0, max(cols)]);
    title('N concentration')
    set(gca,'FontSize',14)
    exportgraphics(gcf,'Project/Figures/cN_it2000_3D.png','ContentType','vector')
    

    set(0, 'currentfigure', 3);
    cla
    view(30,-50);
    cols = state.RN.c;
    plotCellData(model.G, cols);
    colorbar
    caxis([0, max(cols)]);
    title('RN concentration')
    set(gca,'FontSize',14)
    exportgraphics(gcf,'Project/Figures/cRN_it2000_3D.png','ContentType','vector')
    
    drawnow
    pause(0.1);
    
end

%%
% Plot 2d
% Remove empty states (could have been created if solver did not converge)
ind = cellfun(@(state) ~isempty(state), states);
states = states(ind);

figure(1); figure(2); figure(3);

for istate = 2000%1 : numel(states)

    state = states{istate};

    set(0, 'currentfigure', 1);
    cla
    cols = state.R.c;
    plotCellData(model.G, cols);
    colorbar
    caxis([0, max(cols)]);
    title('R concentration')
    set(gca,'FontSize',14)
    exportgraphics(gcf,'Project/Figures/cR_it2000_2D.png','ContentType','vector')
    
    set(0, 'currentfigure', 2);
    cla
    cols = state.N.c;
    plotCellData(model.G, cols);
    colorbar
    caxis([0, max(cols)]);
    title('N concentration')
    set(gca,'FontSize',14)
    exportgraphics(gcf,'Project/Figures/cN_it2000_2D.png','ContentType','vector')

    set(0, 'currentfigure', 3);
    cla
    cols = state.RN.c;
    plotCellData(model.G, cols);
    colorbar
    caxis([0, max(cols)]);
    title('RN concentration')
    set(gca,'FontSize',14)
    exportgraphics(gcf,'Project/Figures/cRN_it2000_2D.png','ContentType','vector')

    drawnow
    pause(0.1);
    
end

%%
% Find out at which point in time the signal can be sent.
for i = 1:n
    if sum(states{i}.R.c.*vols)<sum(states{i}.RN.c.*vols)
        signal = i*dt;
        signal_index = i;
        break
    end
end
% Display result in terminal
signal
signal_index
% Save state at signal to be used to initialize the "clearance" task.
clearance_init_state = states{signal_index};
