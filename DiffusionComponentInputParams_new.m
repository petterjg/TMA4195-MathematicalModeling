classdef DiffusionComponentInputParams_new < ComponentInputParams

    properties
        
        D
        
    end
    
    
    methods

        function paramobj = DiffusionComponentInputParams_new(jsonstruct)
            paramobj = paramobj@ComponentInputParams(jsonstruct);
        end
        
    end
    

    
end
