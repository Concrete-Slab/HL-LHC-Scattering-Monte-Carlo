classdef SecondarySamples
    %SECONDARYSAMPLES Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        secondaries(9,:) double
        secondaryHistories(1,:) ParticleHistory
    end

    properties(Dependent)
        n
    end
    
    methods
        function obj = SecondarySamples(samps,histories)
            %SECONDARYSAMPLES Construct an instance of this class
            arguments
                samps(9,:) double = [];
                histories(1,:) ParticleHistory = ParticleHistory.empty(1,0);
            end
            obj.secondaries = samps;
            obj.secondaryHistories = histories;
        end
        
        function N = get.n(obj)
            N = size(obj.secondaries,2);
        end
    end
end

