classdef ReactionDiffusion_clear < BaseModel

    properties    
        T
        N
        TN
        TNI
        k_on
        k_off
        k_clearance
    end
    
    methods
        
        function model = ReactionDiffusion_clear(paramobj)
            
            model = model@BaseModel();
            
            % All the submodels should have same backend (this is not assigned automaticallly for the moment)
            model.AutoDiffBackend = SparseAutoDiffBackend('useBlocks', false);
            
            %% Setup the model using the input parameters
            fdnames = {'G'            , ...
                       'k_on'            , ...
                       'k_off'            , ...
                       'k_clearance'};
            
            model = dispatchParams(model, paramobj, fdnames);

            model.T = DiffusionComponent_clear(paramobj.T);
            model.N = DiffusionComponent_clear(paramobj.N);
            model.TN = DiffusionComponent_clear(paramobj.TN);
            model.TNI = DiffusionComponent_clear(paramobj.TNI);
            
        end

        
        function model = registerVarAndPropfuncNames_clear(model)
            
            %% Declaration of the Dynamical Variables and Function of the model
            % (setup of varnameList and propertyFunctionList)
            
            model = registerVarAndPropfuncNames_clear@BaseModel(model);
            
            %% Temperature dispatch functions
            fn = @ReactionDiffusion_clear.updateSourceTerm_clear;
            
            inputnames = {{'T', 'c'}, ...
                          {'N', 'c'}, ...
                          {'TN', 'c'}, ...
                          {'TNI', 'c'}};
            model = model.registerPropFunction({{'T', 'source'} , fn, inputnames});
            model = model.registerPropFunction({{'N', 'source'} , fn, inputnames});
            model = model.registerPropFunction({{'TN', 'source'} , fn, inputnames});
            model = model.registerPropFunction({{'TNI', 'source'} , fn, inputnames});
            
        end

        function forces = getValidDrivingForces(model)
            forces = getValidDrivingForces@PhysicalModel(model);
            forces.none = [];
        end
        
        function state = updateSourceTerm_clear(model, state)

            k_on = model.k_on;
            k_off = model.k_off;
            k_clearance = model.k_clearance;
            vols = model.G.cells.volumes;
            
            cT = state.T.c;
            cN = state.N.c;
            cTN = state.TN.c;
            cTNI = state.TNI.c;

            u1 = k_on.*vols.*cT.*cN;
            u2 = k_off.*vols.*cTN;
            u3 = k_clearance.*vols.*cTN;

            state.T.source = -u1 + u2 + u3;
            state.N.source = -u1 + u2;
            state.TN.source = u1 - u2 - u3;
            state.TNI.source = u3;
            
        end
        
        
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            
            state = model.initStateAD(state);

            state.T = model.T.updateFlux_clear(state.T);
            state.N = model.N.updateFlux_clear(state.N);
            state.TN = model.TN.updateFlux_clear(state.TN);
            state.TNI = model.TNI.updateFlux_clear(state.TNI);
            
            state = model.updateSourceTerm_clear(state);

            state.T = model.T.updateMassAccum_clear(state.T, state0.T, dt);
            state.N = model.N.updateMassAccum_clear(state.N, state0.N, dt);
            state.TN = model.TN.updateMassAccum_clear(state.TN, state0.TN, dt);
            state.TNI = model.TNI.updateMassAccum_clear(state.TNI, state0.TNI, dt);
            
            state.T = model.T.updateMassConservation_clear(state.T);
            state.N = model.N.updateMassConservation_clear(state.N);
            state.TN = model.TN.updateMassConservation_clear(state.TN);
            state.TNI = model.TNI.updateMassConservation_clear(state.TNI);
            
            eqs = {}; types = {}; names = {};
            
            eqs{end + 1}   = state.T.massCons;
            names{end + 1} = 'massCons T';
            types{end + 1} = 'cell';
            
            eqs{end + 1}   = state.N.massCons;
            names{end + 1} = 'massCons N';
            types{end + 1} = 'cell';

            eqs{end + 1}   = state.TN.massCons;
            names{end + 1} = 'massCons TN';
            types{end + 1} = 'cell';

            eqs{end + 1}   = state.TNI.massCons;
            names{end + 1} = 'massCons TNI';
            types{end + 1} = 'cell';
                        
            primaryVars = model.getPrimaryVariables();

            %% Setup LinearizedProblem that can be processed by MRST clearton API
            problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
            
        end
        
        
        function state = initStateAD(model, state)
        % initialize a clear cleaned-up state with AD variables

            % initStateAD in BaseModel erase all fields
            clearstate = initStateAD@BaseModel(model, state);
            clearstate.time = state.time;
            state = clearstate;
            
        end 

        
        function primaryvarnames = getPrimaryVariables(model)

            primaryvarnames = {{'T', 'c'}, ...
                               {'N', 'c'}, ...
                               {'TN', 'c'}, ...
                               {'TNI', 'c'}};
            
        end
        

        function model = validateModel(model, varargin)
        % nothing special to do
        end


    end

end
