classdef MCInput
    %MCIN Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        samples(9,:) double
        geometry(1,1) Geometry = OffsetRectangle;
        material(1,1) Material = Material
    end

    properties(Dependent,SetAccess=private)
        nSamples(1,1) double
    end
    
    methods % constructor
        function obj = MCInput(geometry,material,options)
            %MCIN Construct an instance of this class
            %   Detailed explanation goes here
            arguments
                geometry(1,1) Geometry = OffsetRectangle
                material(1,1) Material = Material
                options.n(1,1) double {mustBeInteger,mustBePositive}
                options.samples(9,:) double 
                options.echo(1,1) string {mustBeMember(options.echo,["off","on"])} = "off";
            end
            if (isfield(options,'samples'))
                obj.samples = options.samples;
            else
                if ~isfield(options,'n')
                    n = 0;
                else
                    n = options.n;
                end
%                 warning("Generating %d samples, expect performance loss",n);
                samps = geometry.generate(n,"echo",options.echo);
                obj.samples = samps;
            end
            obj.geometry = geometry;
            obj.material = material;
            
        end
    end

    methods
        function view(obj)
%             if nargin == 1
%                 f = figure;
%             end
            % fix integration with figure input




            obj.geometry.view
            hold on
            samps = obj.samples(1:3,:);
            plot3(samps(1,:),samps(2,:),samps(3,:),"r.")
            theta = linspace(-pi,pi,250);
            circVals = Consts.TubeRadius*[cos(theta);sin(theta)];
            %plot(circVals(1,:),circVals(2,:),"k")
            hold off
        end
    end

    methods % getters
        function n = get.nSamples(obj)
            n = size(obj.samples,2);
        end
    end
end

