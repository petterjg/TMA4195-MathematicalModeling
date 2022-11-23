classdef ReactionDiffusion_new < BaseModel

    properties    
        R
        N
        RN
        k_on
        k_off
    end
    
    methods
        
        function model = ReactionDiffusion_new(paramobj)
            
            model = model@BaseModel();
            
            % All the submodels should have same backend (this is not assigned automaticallly for the moment)
            model.AutoDiffBackend = SparseAutoDiffBackend('useBlocks', false);
            
            %% Setup the model using the input parameters
            fdnames = {'G'            , ...
                       'k_on'            , ...
                       'k_off'};
            
            model = dispatchParams(model, paramobj, fdnames);

            model.R = DiffusionComponent_new(paramobj.R);
            model.N = DiffusionComponent_new(paramobj.N);
            model.RN = DiffusionComponent_new(paramobj.RN);
            
        end

        
        function model = registerVarAndPropfuncNames_new(model)
            
            %% Declaration of the Dynamical Variables and Function of the model
            % (setup of varnameList and propertyFunctionList)
            
            model = registerVarAndPropfuncNames_new@BaseModel(model);
            
            %% Temperature dispatch functions
            fn = @ReactionDiffusion_new.updateSourceTerm_new;
            
            inputnames = {{'R', 'c'}, ...
                          {'N', 'c'}};
            model = model.registerPropFunction({{'R', 'source'} , fn, inputnames});
            model = model.registerPropFunction({{'N', 'source'} , fn, inputnames});
            model = model.registerPropFunction({{'RN', 'source'} , fn, inputnames});
            
        end

        function forces = getValidDrivingForces(model)
            forces = getValidDrivingForces@PhysicalModel(model);
            forces.none = [];
        end
        
        function state = updateSourceTerm_new(model, state)

            k_on = model.k_on;
            k_off = model.k_off;
            vols = model.G.cells.volumes;
            
            cR = state.R.c;
            cN = state.N.c;
            cRN = state.RN.c;

            u1 = k_on.*vols.*cR.*cN;
            u2 = k_off.*vols.*cRN;

            state.R.source = -u1 + u2;
            state.N.source = -u1 + u2;
            state.RN.source = u1 - u2;
            
        end
        
        
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            
            state = model.initStateAD(state);

            state.R = model.R.updateFlux_new(state.R);
            state.N = model.N.updateFlux_new(state.N);
            state.RN = model.RN.updateFlux_new(state.RN);
            
            state = model.updateSourceTerm_new(state);

            state.R = model.R.updateMassAccum_new(state.R, state0.R, dt);
            state.N = model.N.updateMassAccum_new(state.N, state0.N, dt);
            state.RN = model.RN.updateMassAccum_new(state.RN, state0.RN, dt);
            
            state.R = model.R.updateMassConservation_new(state.R);
            state.N = model.N.updateMassConservation_new(state.N);
            state.RN = model.RN.updateMassConservation_new(state.RN);
            
            eqs = {}; types = {}; names = {};
            
            eqs{end + 1}   = state.R.massCons;
            names{end + 1} = 'massCons R';
            types{end + 1} = 'cell';
            
            eqs{end + 1}   = state.N.massCons;
            names{end + 1} = 'massCons N';
            types{end + 1} = 'cell';

            eqs{end + 1}   = state.RN.massCons;
            names{end + 1} = 'massCons RN';
            types{end + 1} = 'cell';
                        
            primaryVars = model.getPrimaryVariables();

            %% Setup LinearizedProblem that can be processed by MRST Newton API
            problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
            
        end
        
        
        function state = initStateAD(model, state)
        % initialize a new cleaned-up state with AD variables

            % initStateAD in BaseModel erase all fields
            newstate = initStateAD@BaseModel(model, state);
            newstate.time = state.time;
            state = newstate;
            
        end 

        
        function primaryvarnames = getPrimaryVariables(model)

            primaryvarnames = {{'R', 'c'}, ...
                               {'N', 'c'}, ...
                               {'RN', 'c'}};
            
        end
        

        function model = validateModel(model, varargin)
        % nothing special to do
        end


    end

end
