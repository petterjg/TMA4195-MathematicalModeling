classdef ReactionDiffusionInputParams_clear < InputParams

    properties
        G      
        T
        N
        TN
        TNI
        k_on
        k_off
        k_clearance
    end
    
    methods
        
        function paramobj = ReactionDiffusionInputParams_clear(jsonstruct)
            
            paramobj = paramobj@InputParams(jsonstruct);
            
            pick = @(fd) pickField(jsonstruct, fd);
            
            paramobj.T = DiffusionComponentInputParams_clear(pick('T'));
            paramobj.N = DiffusionComponentInputParams_clear(pick('N'));
            paramobj.TN = DiffusionComponentInputParams_clear(pick('TN'));
            paramobj.TNI = DiffusionComponentInputParams_clear(pick('TNI'));

            paramobj = paramobj.validateInputParams_clear();
            
        end

        function paramobj = validateInputParams_clear(paramobj)

            if ~isempty(paramobj.G)
                paramobj.T.G = paramobj.G;
                paramobj.N.G = paramobj.G;
                paramobj.TN.G = paramobj.G;
                paramobj.TNI.G = paramobj.G;
            end
            
        end
        
    end
    
end
