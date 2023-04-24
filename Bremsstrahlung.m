classdef Bremsstrahlung < SoftProcess
    %BREMSSTRAHLUNG Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        X0
    end
    
    methods
        function obj = Bremsstrahlung(material,particle)
            %BREMSSTRAHLUNG Construct an instance of this class
            %   Detailed explanation goes here
            arguments
                material(1,1) Material
                particle(1,:) ParticleHandle {mustBeScalarOrEmpty} = ParticleHandle.empty
            end
            obj.X0 = material.X0;
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.X0 + inputArg;
        end
    end
end

