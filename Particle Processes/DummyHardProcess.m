classdef DummyHardProcess < HardProcess
    %DUMMYHARDPROCESS Dummy hard process that has no probability of
    % occuring and does not affect the particle
    %   Instantiated when creating an empty array HardProcess array
    
    properties(Dependent) % inherited
        sigmaTot
    end
    
    methods
        function obj = DummyHardProcess(material)
            %DUMMYHARDPROCESS Construct an instance of this class
            %   Detailed explanation goes here
            arguments
                material(1,1) Material = Material;
            end
            obj@HardProcess(material)
        end
        
        function st = get.sigmaTot(obj)
            st = 0;
        end

        function [delE,delTheta,secondaryParticle] = interact(obj)
            delE = 0;
            delTheta = [0;1;0;1];
            secondaryParticle = double.empty(6,0);
        end
    end
end

