classdef DummySecondaryProcess < HardProcess
    %DUMMYSECONDARYPROCESS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Dependent)
        sigmaTot
    end
    
    methods
        function obj = DummySecondaryProcess(material)
            %DUMMYSECONDARYPROCESS Construct an instance of this class
            arguments
                material(1,1) Material
            end
            obj@HardProcess(material)
        end
        
        function st = get.sigmaTot(obj)
            st = 10/obj.n;
            % path length of 1m results from this value
        end

        function [dE,dTheta,secondaryParticle] = interact(obj)
            dE = 0;
            z = rand(2,1);
            dTheta = [sin(z(1));cos(z(1));sin(0.05*z(2));cos(0.05*z(2))];
%             dTheta = 0;
            secenergy = obj.particle.energy;
            seccharge = obj.particle.charge;
            dir = [2;2;-0.5].*rand(3,1)-[1;1;-1];
%             dir = rand(3,1);
            dir = dir./norm(dir);
            secmom = obj.particle.momentum*dir;
            secondaryParticle = [obj.particle.energy;secmom;obj.particle.energy;obj.particle.mass];
        end
    end
end

