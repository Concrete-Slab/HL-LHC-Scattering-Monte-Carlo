classdef(Abstract) HardProcess < matlab.mixin.Heterogeneous
    %UNTITLED Summary of this class goes here
    %   Assume isotropic material...

    properties(Dependent,SetAccess=private)
        lambda(1,1) double {mustBeNonnegative}      % mean free path, m
        M0(1,1) double {mustBeNonnegative}          % inverse mean free path, /m
    end

    properties(Dependent,Abstract)  % must be implemented
        sigmaTot(1,1) double {mustBeNonnegative}    % total interaction cross section, m^2
    end

    properties % final
        n(1,1) double
        particle(1,:) ParticleHandle {mustBeScalarOrEmpty}
    end

    methods
        function obj = HardProcess(material,particlehandle)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            arguments
                material(1,1) Material
                particlehandle(1,:) ParticleHandle {mustBeScalarOrEmpty} = ParticleHandle.empty(1,0)
            end
            obj.n = material.n;
            obj.particle = particlehandle;
        end
    end

    methods % getters

        function l = get.lambda(obj)
            l = 1./(obj.sigmaTot*obj.n);
        end

        function il = get.M0(obj)
            il = obj.sigmaTot*obj.n;
        end
    end

    methods(Abstract) % abstract inverse random sampling functions
        [delE,delTheta,secondaryParticle] = interact(obj)
    end

    methods(Static,Access=protected) % overload
        function hp = getDefaultScalarElement()
            hp = DummyHardProcess;
        end
    end
end

