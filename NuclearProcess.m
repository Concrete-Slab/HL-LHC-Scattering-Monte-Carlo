classdef NuclearProcess < HardProcess
    %NuclearProcess Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        A(1,1) double       % average relative atomic mass
    end

    properties(Dependent)
        sigmaTot
        sigmappElastic
        sigmappSD
        
    end
    
    methods
        function obj = NuclearProcess(material,particle)
            arguments
                material(1,1) Material = Material
                particle(1,:) ParticleHandle {mustBeScalarOrEmpty} = ParticleHandle.empty
            end
            %NuclearProcess Construct an instance of this class
            %   Detailed explanation goes here
            obj@HardProcess(material,particle);
            obj.A = sum(material.wt.*material.A);
        end
        
        function st = get.sigmaTot(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            
        end

        function [dE,dAngle,newSecondaryInfo] = interact(obj)
            newSecondaryInfo = double.empty(6,0);
            
        end
    end

    methods
        
    end
end

