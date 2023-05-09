classdef Bremsstrahlung < SoftProcess
    %BREMSSTRAHLUNG
    %   Continuous radiative energy loss by electrons
    
    properties(SetAccess=private)
        X0
    end

    properties(Dependent,SetAccess=protected)
        dEdx
    end
    
    methods
        function obj = Bremsstrahlung(material,particle)
            %BREMSSTRAHLUNG Construct an instance of this class
            %   Detailed explanation goes here
            arguments
                material(1,1) Material = "Chromox"
                particle(1,:) ParticleHandle {mustBeScalarOrEmpty} = ParticleHandle.empty
            end
            obj.X0 = material.X0;
            obj.particle = particle;
        end
        
        function [dE,dAngles,dPos] = update(obj,delX)
            %UPDATE Apply an energy change due to brehmsstrahlung
            % only apply update to electrons
            if obj.particle.charge == -1 && obj.particle.mass == Consts.megev
                dE = obj.dEdx * delX;
            else
                dE = 0;
            end
            dAngles = 0;
            dPos = 0;
        end
    end

    methods % getters
        function e = get.dEdx(obj)
            e = -1 * obj.particle.energy/obj.X0;
        end
    end
end

