classdef ParticleHandle < handle
    %PARTICLEHANDLE Stores particle properties as a handle
    % The handle class allows changes to the particle to be reflected
    % across all models, enabling updates from one process to affect all
    % other processes
    
    properties
        energy      % Total particle energy in GeV
    end

    properties(SetAccess=private)
        charge      % Particle charge in e
        mass        % Particle mass in natural units, GeV/c^2
    end

    properties(Dependent)
        momentum    % Particle momentum in natural units, GeV/c
    end

    properties(Dependent,SetAccess=private)
        gamma       % gamma factor
        beta        % beta factor - fraction of c
        kineticEnergy
    end
    
    methods % constructor
        function obj = ParticleHandle(energy,charge,mass)
            %PARTICLEHANDLE Constructor
            arguments
                energy(1,1) double      % GeV
                charge(1,1) double      % e
                mass(1,1) double        % GeV/c^2
            end
            obj.energy = energy;
            obj.charge = charge;
            obj.mass = mass;
        end
    end

    methods % getters/setters
        function p = get.momentum(obj)
            p = sqrt(obj.energy^2-obj.mass^2);
        end

        function set.momentum(obj,p)
            obj.energy = sqrt(obj.mass^2+p^2);
        end

        function y = get.gamma(obj)
            y = obj.energy/obj.mass;
        end

        function b = get.beta(obj)
            b = sqrt(1-1/obj.gamma^2);
        end

        function ek = get.kineticEnergy(obj)
            ek = obj.mass*(obj.gamma-1);
        end
    end
end

