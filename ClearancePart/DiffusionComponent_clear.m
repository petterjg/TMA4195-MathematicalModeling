classdef DiffusionComponent_clear < BaseModel

    properties
        
        D % Diffusion coefficient
    end

    methods
        
        function model = DiffusionComponent_clear(paramobj)

            model = model@BaseModel();
            
            % OBS : All the models should have same backend (this is not assigned automaticallly for the moment)
            model.AutoDiffBackend = SparseAutoDiffBackend('useBlocks', false);
            
            fdnames = {'G', ...
                       'D'};
            
            model = dispatchParams(model, paramobj, fdnames);

            model.operators = localSetupOperators(model.G);
            

        end

        function model = registerVarAndPropfuncNames_clear(model)
        %% Declaration of the Dynamical Variables and Function of the model
        %  Note : this is more like a documentation and is not used in assembly

            model = registerVarAndPropfuncNames_clear@BaseModel(model);
            
            varnames = {'c'       , ...
                        'accum'   , ...
                        'source'  , ...
                        'flux'    , ...
                        'massCons'};
            
            model = model.registerVarNames(varnames);

            fn = @DiffusionComponent_clear.updateFlux_clear;
            inputnames = {'c'};
            model = model.registerPropFunction({'flux', fn, inputnames});

            fn = @DiffusionComponent_clear.updateMassConservation_clear;
            inputnames = {'accum', 'flux', 'source'};
            model = model.registerPropFunction({'massCons', fn, inputnames});
            
            fn = @DiffusionComponent_clear.updateMassAccum_clear;
            inputnames = {'c'};
            model = model.registerPropFunction({'accum', fn, inputnames});

        end

        
        function state = updateFlux_clear(model, state)
        % Assemble electrical current which is stored in :code:`state.j` 

            D  = model.D;
            op = model.operators;
            T  = op.T;
            
            c   = state.c;
            
            flux = - D.*T.*op.Grad(c);

            state.flux = flux;
            
        end

        function state = updateMassAccum_clear(model, state, state0, dt)

            vols = model.G.cells.volumes;
            
            c  = state.c;
            c0 = state0.c;

            state.accum = vols.*(c - c0)/dt;
            
        end
        
        
        function state = updateMassConservation_clear(model, state)
        % Assemble residual of the charge conservation equation which is stored in :code:`state.chargeCons`

            accum    = state.accum;
            flux     = state.flux;
            source   = state.source;
            
            massCons = assembleConservationEquation(model, flux, 0, source, accum);
            
            state.massCons = massCons;
            
        end
        
    end
    
end


