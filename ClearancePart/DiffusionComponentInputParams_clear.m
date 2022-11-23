classdef DiffusionComponentInputParams_clear < ComponentInputParams

    properties
        
        D
        
    end
    
    
    methods

        function paramobj = DiffusionComponentInputParams_clear(jsonstruct)
            paramobj = paramobj@ComponentInputParams(jsonstruct);
        end
        
    end
    

    
end
