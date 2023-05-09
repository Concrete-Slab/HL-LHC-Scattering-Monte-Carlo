classdef SymRectangles < Geometry
    %HORIZONTALSYMMETRICGAP Summary of this class goes here
    %   Detailed explanation goes here

    properties(SetAccess=private)
        width(1,1) double {mustBePositive} = 0.01
        height(1,1) double {mustBePositive} = 0.01
        centralGap(1,1) double {mustBeNonnegative} = 0
    end

    properties(Dependent,SetAccess=private)
        xminus
        xplus
    end

    methods % geometrical property getters
        function xminus = get.xminus(obj)
            xminus = max(obj.xv(1:5));
        end

        function xplus = get.xplus(obj)
            xplus = min(obj.xv(7:11));
        end
    end

    methods % constructor
        function obj = SymRectangles(width,height,centralGap,thickness)
            arguments
                width(1,1) double {mustBePositive} = 0.01
                height(1,1) double {mustBePositive} = 0.01
                centralGap(1,1) double {mustBeNonnegative} = 0.001
                thickness(1,1) double {mustBePositive} = 0.001
            end

            %HORIZONTALSYMMETRICGAP Construct an instance of this class
            %   Detailed explanation goes here
            leftRectX = [(-centralGap/2), (-centralGap/2-width), (-centralGap/2-width), (-centralGap/2), (-centralGap/2)];
            leftRectY = [-height/2, -height/2, height/2, height/2, -height/2];
            rightRectX = [(centralGap/2 + width), (centralGap/2), (centralGap/2), (centralGap/2 + width), (centralGap/2 + width)];
            rightRectY = leftRectY;
            xv = [leftRectX,NaN,rightRectX];
            yv = [leftRectY,NaN,rightRectY];
            obj@Geometry(xv,yv,0,thickness);
            obj.width = width;
            obj.height = height;
            obj.centralGap = centralGap;
        end
    end


%     methods % abstract inherited
%         function b = isInside(obj,r)
%             arguments(Input)
%                 obj(1,1) Geometry
%                 r(3,:) double
%             end
%             arguments(Output)
%                 b(1,:) logical
%             end
%             x = r(1,:);
%             y = r(2,:);
%             z = r(3,:);
%             % draw up conditions that must be met for inside
%             xConds = or(and(x<=obj.xmax,x>=obj.xplus),and(x>=obj.xmin,x<=obj.xminus));
%             yConds = and(y<=obj.ymax,y>=obj.ymin);
%             zConds = and(z<=obj.zmax(y),z>=obj.zmin(y));
% 
%             b = and(and(xConds,yConds),zConds);
%         end
% 
%         function r = edgeIntersect(obj,pos,direction)
%             arguments(Input)
%                 obj(1,1) Geometry
%                 pos(3,:) double {mustBeFinite,Geometry.mustBeInside(pos,obj)}
%                 direction(3,:) double {mustBeFinite,Geometry.mustBeSameLength(direction,pos)}
%             end
%             arguments(Output)
%                 r(3,:) double {mustBeFinite}
%             end
%             % define number of input sample positions
%             N = size(pos,2);
% 
%             % input validation
%             b = obj.isInside(pos);
%             if ~b
%                 error("Sample was not inside volume when edgeIntersect called")
%             end
% 
%             directionMags = sqrt(diag(direction'*direction))';
% 
%             direction = direction./directionMags;
% 
%             % forwards z
%             ratios = zeros(8,N);
%             ratios(1,:) = (obj.xmax-pos(1,:))./direction(1,:);
%             ratios(2,:) = (obj.xplus-pos(1,:))./direction(1,:);
%             ratios(3,:) = (obj.xminus-pos(1,:))./direction(1,:);
%             ratios(4,:) = (obj.xmin-pos(1,:))./direction(1,:);
%             ratios(5,:) = (obj.ymax-pos(2,:))./direction(2,:);
%             ratios(6,:) = (obj.ymin-pos(2,:))./direction(2,:);
%             ratios(7,:) = (obj.zmin(pos(2,:))-pos(3,:))./direction(3,:);
%             ratios(8,:) = (obj.zmax(pos(2,:))-pos(3,:))./direction(3,:);
%             minRatios = zeros(1,N);
%             for i=1:N
%                 currRatios = ratios(:,i);
%                 currRatios = currRatios(abs(currRatios)<Inf);
%                 I = find(currRatios==0);
%                 if ~isempty(I)
%                     consideredVals = zeros(length(I));
%                     for j = 1:length(I)
%                         if mod(I(j),2) == 0
%                             otherIndex = I(j)-1;
%                         else
%                             otherIndex = I(j)+1;
%                         end
%                         if currRatios(otherIndex)>0
%                             consideredVals(j) = currRatios(otherIndex);
%                         else
%                             consideredVals(j) = 0;
%                         end   
%                     end
%                     currRatios(I) = consideredVals;
%                     currentMin = min(currRatios(currRatios>=0));
%                 else
%                     currentMin = min(currRatios(currRatios>0));
%                 end
% 
%                 if isempty(currentMin) % if sample was originally not inside the volume
%                     % do nothing might not be the right course of action...
%                     currentMin = 0;
%                 end
% 
%                
%                 minRatios(i) = currentMin;
%             end
% 
%             % FIX COMPUTATIONAL ROUNDING ERROR
%             dp = floor(min(min(-log10(eps(pos)))));
%             
%             r = round(pos + direction.*minRatios,dp);
%         end
% 
%         % TODO not functional - might be functional now
% 
%         function [b,reEntryPos] = willReEnter(obj,pos,direction)
%             arguments(Input)
%                 obj(1,1) SymRectangles
%                 pos(3,:) double
%                 direction(3,:) double
%             end
%             arguments(Output)
%                 b(1,:) logical
%                 reEntryPos(3,:) double
%             end
%             
%             % only way a particle can re-enter is through xplus or xminus
%             useplus = and(direction(1,:)>0,pos(1,:)<obj.xplus);
%             useminus = and(direction(1,:)<0,pos(1,:)>obj.xminus);
%             dp = floor(min(min(-log10(eps(pos)))));
%             b = false(1,size(pos,2));
%             reEntryPlus = round(pos(:,useplus) + (obj.xplus-pos(1,useplus))./direction(1,useplus).*direction(:,useplus),dp);
%             b(useplus) = obj.isInside(reEntryPlus);
%             reEntryMinus = round(pos(:,useminus) + (obj.xminus-pos(1,useminus))./direction(1,useminus).*direction(:,useminus),dp);
%             b(useminus) = obj.isInside(reEntryMinus);
%             reEntryPos = zeros(3,length(b));
%             reEntryPos(:,useplus) = reEntryPlus;
%             reEntryPos(:,useminus) = reEntryMinus;
%             reEntryPos = reEntryPos(:,b);
%         end
% 
%     end
end


