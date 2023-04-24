classdef DummySoftProcess < SoftProcess
    %DUMMYSOFTPROCESS Summary of this class goes here
    %   Detailed explanation goes here

    properties(Dependent,SetAccess=protected)
        dEdx
    end
    
    methods
        function obj = DummySoftProcess
            %DUMMYSOFTPROCESS Construct an instance of this class
            %   Detailed explanation goes here
        end

        function e = get.dEdx(obj)
            e = 0;
        end
        
        function [delE,delTheta,delPos] = update(obj,delX)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            arguments(Output)
                delE(1,1) double
                delTheta(2,1) double
                delPos(2,1) double
            end
            delE = 0;
            delTheta = zeros(2,1);
            delPos = zeros(2,1);
        end
    end
end

