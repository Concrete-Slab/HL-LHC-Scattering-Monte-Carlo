classdef(Abstract) SoftProcess < matlab.mixin.Heterogeneous
    %SOFTPROCESS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        particle(1,:) ParticleHandle {mustBeScalarOrEmpty} = ParticleHandle.empty
    end

    properties(Dependent,SetAccess=protected,Abstract)
        dEdx(1,1) double
    end

    
    methods
        function obj = SoftProcess(particlehandle)
            %SOFTPROCESS Construct an instance of this class
            %   Detailed explanation goes here
            arguments
                particlehandle(1,:) ParticleHandle {mustBeScalarOrEmpty} = ParticleHandle.empty
            end
            obj.particle = particlehandle;
        end
    end

    methods(Abstract)
        [dE,dAngles,dPos] = update(obj,dr);
%         dx = minStepLength(obj,initialConditions)
    end

    methods(Static,Access=protected) % overload
        function sp = getDefaultScalarElement()
            sp = DummySoftProcess;
        end
    end
    
end

