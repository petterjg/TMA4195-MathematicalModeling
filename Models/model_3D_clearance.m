close all

mrstModule add ad-core mrst-gui 

jsonfile = fileread('input_clear.json');
jsonstruct = jsondecode(jsonfile);

paramobj = ReactionDiffusionInputParams_clear(jsonstruct);


% Define grid
radius = 0.22*micro*meter;
height = 15*nano*meter;
radiusLayers = 10;
heightLayers = 10;

% Initialize variables
avo = 6.022e23;
% MUST RUN AFTER model_3D_and_2D with the same dimension has run.
dim = 3;
if dim == 3
    % 3D
    Grid = CylinderGrid(radius,height,radiusLayers,heightLayers);
    vols = Grid.cells.volumes;
    nc = Grid.cells.num;
    topLayerIndices = (1 : nc/heightLayers);
    bottomLayerIndices = (nc - nc/heightLayers + 1 : nc);
    % same as receptors? Should be enough glia cells to react with a
    % considerable share of Ns.
    cT0 = 1e6/(avo*micro^2)/(radius/radiusLayers);
end
if dim == 2
    % 2D
    Grid = CircleGrid(radius,radiusLayers);
    vols = Grid.cells.volumes*height;
    nc = Grid.cells.num;
    topLayerIndices = (1 : nc);
    bottomLayerIndices = topLayerIndices;
    cT0 = 1000/(avo*micro^2)/(height);
end

% Get initial concentration of N from point in time the signal was sent
% after running model_3D.m.
% Add concentration of RN to N as all Ns are released on signal.
cN0 = clearance_init_state.N.c + clearance_init_state.RN.c;

paramobj.k_on = 1e5;
paramobj.k_off = 5;
paramobj.k_clearance = 50;
paramobj.N.D = 3e-10;
paramobj.T.D = 0;
paramobj.TN.D = 0;
paramobj.TNI.D = 3e-10;
paramobj.G = Grid;

paramobj = paramobj.validateInputParams_clear();

model = ReactionDiffusion_clear(paramobj);


% Setup schedule
total = 1e-3;
n  = 1000;
dt = total/n;
step = struct('val', dt*ones(n, 1), 'control', ones(n, 1));

control.none = [];
schedule = struct('control', control, 'step', step);

% Setup initial state
% Density of recptors: Density on membrane * volume of cells on bottom
% level / height on cells on bottom level


initcase = 1;
switch initcase
  case 1
    cT = zeros(nc, 1);
    % Place T only in the outermost cylinderlayer
    cT(find((vols > vols(20*radiusLayers)*0.99).*(vols < vols(20*radiusLayers)*1.01)==1)) = cT0;
    cN = zeros(nc, 1);
    %cN(topLayerIndices) = cN0;
    cN = cN0;
    cTN = zeros(nc, 1);
    cTNI = zeros(nc, 1);
  case 2
    cT = ones(nc, 1);
    cN = ones(nc, 1);
    cTN = zeros(nc, 1);
    cTNI = zeors(nc, 1);
end

initstate.T.c = cT;
initstate.N.c = cN;
initstate.TN.c = cTN;
initstate.TNI.c = cTNI;

% run simulation

nls = NonLinearSolver();
nls.errorOnFailure = false;

[~, states, report] = simulateScheduleAD(initstate, model, schedule, 'NonLinearSolver', nls);

%%

% Remove empty states (could have been created if solver did not converge)
ind = cellfun(@(state) ~isempty(state), states);
states = states(ind);

figure(1); figure(2); figure(3); figure(4);

for istate = 610%1 : numel(states)

    state = states{istate};

    set(0, 'currentfigure', 1);
    cla
    view(30,-50);
    cols = state.T.c;
    plotCellData(model.G, cols);
    colorbar
    caxis([0, max(cols)]);
    title('T concentration')
    set(gca,'FontSize',14)
    exportgraphics(gcf,'Project/Figures/cTc_it610_3D.png','ContentType','vector')
    
    set(0, 'currentfigure', 2);
    cla
    view(30,-50);
    cols = state.N.c;
    plotCellData(model.G, cols);
    colorbar
    caxis([0, max(cols)]);
    title('N concentration')
    set(gca,'FontSize',14)
    exportgraphics(gcf,'Project/Figures/cNc_it610_3D.png','ContentType','vector')

    set(0, 'currentfigure', 3);
    cla
    view(30,-50);
    cols = state.TN.c;
    plotCellData(model.G, cols);
    colorbar
    caxis([0, max(cols)]);
    title('TN concentration')
    set(gca,'FontSize',14)
    exportgraphics(gcf,'Project/Figures/cTNc_it610_3D.png','ContentType','vector')

    set(0, 'currentfigure', 4);
    cla
    view(30,-50);
    cols = state.TNI.c;
    plotCellData(model.G, cols);
    colorbar
    caxis([0, max(cols)]);
    title('N_{inactive} concentration')
    set(gca,'FontSize',14)
    exportgraphics(gcf,'Project/Figures/cTNIc_it610_3D.png','ContentType','vector')

    drawnow
    pause(0.1);
    
end

%%
% Find out at which point in time the synapse is cleared.
proportion = 0.99; % How big share of the neurotransmitters must have become inactive
for i = 1:n
    if sum(states{i}.TNI.c.*vols)>proportion*(sum(states{i}.TNI.c.*vols) + sum(states{i}.N.c.*vols))
        clearance = i*dt;
        clearance_index = i;
        break
    end
end
% Display result in terminal
clearance
clearance_index