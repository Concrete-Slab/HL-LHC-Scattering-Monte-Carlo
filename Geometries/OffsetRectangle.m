classdef OffsetRectangle < Geometry
    %OFFSETRECTANGLE Summary of this class goes here
    %   Detailed explanation goes here

    properties
        width
        height
        offset
    end

    properties(Constant)
        defaultWidth = 0.01;
        defaultHeight = 0.01;
        defaultOffset = Consts.MaxSigma*4.8+Consts.MaxMu;
    end

    methods
        function obj = OffsetRectangle(width,height,offset,thickness,xRotation)
            %OFFSETRECTANGLE Construct an instance of this class
            %   Detailed explanation goes here
            arguments
                width(1,1) double {mustBePositive} = 0.01
                height(1,1) double {mustBePositive} = 0.01
                offset(1,1) double {mustBePositive} = 0.01;
                thickness(1,1) double {mustBePositive} = 0.001;
                xRotation(1,1) double {mustBeNonnegative,mustBeLessThan(xRotation,90)} = 0;
            end
            xv = [offset+width offset offset offset+width offset+width];
            yv = [-height/2 -height/2 height/2 height/2 -height/2];
            obj@Geometry(xv,yv,xRotation,thickness)
            obj.width = width;
            obj.height = height;
            obj.offset = offset;
        end
    end
end