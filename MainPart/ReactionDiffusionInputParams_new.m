classdef ReactionDiffusionInputParams_new < InputParams

    properties
        G      
        R
        N
        RN
        k_on
        k_off
    end
    
    methods
        
        function paramobj = ReactionDiffusionInputParams_new(jsonstruct)
            
            paramobj = paramobj@InputParams(jsonstruct);
            
            pick = @(fd) pickField(jsonstruct, fd);
            
            paramobj.R = DiffusionComponentInputParams_new(pick('R'));
            paramobj.N = DiffusionComponentInputParams_new(pick('N'));
            paramobj.RN = DiffusionComponentInputParams_new(pick('RN'));

            paramobj = paramobj.validateInputParams_new();
            
        end

        function paramobj = validateInputParams_new(paramobj)

            if ~isempty(paramobj.G)
                paramobj.R.G = paramobj.G;
                paramobj.N.G = paramobj.G;
                paramobj.RN.G = paramobj.G;
            end
            
        end
        
    end
    
end
