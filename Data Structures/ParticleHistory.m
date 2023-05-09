classdef ParticleHistory
    %PARTICLEHISTORY Summary of this class goes here
    %   Detailed explanation goes here
    

    properties(SetAccess=private)
        positions(3,:) double {mustBeFinite}
    end

    properties(Access=private)
        defaultList(3,:) double
        increment(1,1) double
    end

    properties(SetAccess=private)
        n(1,1) double {mustBeInteger,mustBeNonnegative}
    end
    
    methods
        function obj = ParticleHistory(estimatedLength,initialPosition)
            %PARTICLEHISTORY Construct an instance of this class
            %   Detailed explanation goes here
            arguments
                estimatedLength(1,1) double {mustBeInteger,mustBePositive} = 3
                initialPosition(3,:) double = double.empty(3,0);
            end
            obj.defaultList = zeros(3,estimatedLength);
            obj.positions = obj.defaultList;
            if isempty(initialPosition)
                obj.increment = 0;
            else
                obj.increment = 1;
                obj.positions(:,1) = initialPosition;
            end
        end

        function obj = write(obj,newPos)
            arguments
                obj(1,1) ParticleHistory
                newPos(3,:) double {mustBeFinite}
            end
            oldIncrement = obj.increment + 1;
            obj.increment = obj.increment + size(newPos,2);
            if obj.increment > size(obj.positions,2)
                newAmount = ceil((obj.increment-size(obj.positions,2))./size(obj.defaultList,2));
                obj.positions = [obj.positions, repmat(obj.defaultList,1,newAmount)];
            end
            obj.positions(:,oldIncrement:obj.increment) = newPos;
        end

        function obj = memsave(obj)
            obj.positions = obj.positions(:,1:obj.increment);
        end
        
        function n = get.n(obj)
            n = size(obj.positions,2);
        end
    end
end

