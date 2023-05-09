classdef Geometry < handle & matlab.mixin.Heterogeneous
    %GEOMETRY Summary of this class goes here
    %   Detailed explanation goes here

    properties(SetAccess=private)                       % all lengths in metres
        xv(1,:) double                                  % x positions of vertices
        yv(1,:) double                                  % y positions of vertices
        integratedProbability(1,:) double {mustBeNonnegative,mustBeScalarOrEmpty}
        zvMin(1,:) double                               % z positions of upstream vertices
        zvMax(1,:) double                               % z positions of downstream vertices
        area(1,1) double                                % frontal area of screen
        precision(1,1) double {mustBeInteger,mustBePositive} = 16 % maximum number of decimal places to round to
    end

    properties % public settable
        pdf(:,1) function_handle {mustBeScalarOrEmpty}  % proton pdf, equal in x and y dir
        xRotation(1,1) double {mustBeFinite}            % angle between surface normal and beam direction, degrees
        thickness(1,1) double {mustBePositive} = 0.001  % screen thickness, m
    end

    properties(Dependent,SetAccess=protected) % xv/yv dependent
        xmax(1,1) double
        xmin(1,1) double
        ymax(1,1) double
        ymin(1,1) double
        xLimits(2,:) double % row 1: lower x limits for each volume, row 2: upper x limits for each volume
        yLimits(2,:) double % similar structure to xLimits
    end

    properties(Dependent,SetAccess=protected) % cdf for proton position in x and y directions
        % cdf properties are dependent on pdf and the runtime cdfResolution
        % requires specific implementation of a getter by the subclass
        % row 1 lists tabulated position
        % row 2 lists tabulated probability P(X<x)
        cdfx(2,:) double {mustBeNonnegative}        % tabulated cdf in x dir
        cdfy(2,:) double {mustBeNonnegative}        % tabulated cdf in y dir   
    end

    methods(Access=protected) % constructor
        function obj = Geometry(xv,yv,xRotation,thickness)
            %GEOMETRY Construct an instance of this class
            %   Detailed explanation goes here
            arguments
                xv(1,:) double {mustBeNonempty}
                yv(1,:) double {mustBeNonempty}
                xRotation(1,1) double {mustBeFinite} = 0
                thickness(1,1) double {mustBePositive} = 0.001;
            end
            r = sqrt(xv.^2 + yv.^2);
            if any(r>=Consts.TubeRadius)
                warning("Some vertices of the screen are outside the beam tube")
            end
            obj.xv = xv;
            obj.yv = yv;
            obj.thickness = thickness;
            obj.xRotation = xRotation;
            

            % calculate frontal area from vertices - shoelace formula
            nVertices = length(xv)-1;
            obj.area = 1/2 * abs(sum((yv(1:nVertices)+yv(2:nVertices+1)).*(xv(1:nVertices)-xv(2:nVertices+1))));

            
        end
    end

    methods % public, simulation
        function samples = generate(obj,n,options)
            % SAMPLES
            % Samples are stored as column vectors of length 8
            % entries 1-3 store 3d position
            % entries 4-7 store 3d four momentum  in natural units GeV/c
            % entry 8 stores particle charge
            % entry 9 stores rest mass in natural units GeV/c^2 for numerical accuracy
            arguments(Input)
                obj(1,1) Geometry
                n(1,1) double {mustBeNonnegative,mustBeInteger}
                options.type(1,1) char {mustBeMember(options.type,{'p','e'})} = 'p'
                options.echo(1,1) string {mustBeMember(options.echo,["on","off"])} = "on"
            end
            arguments(Output)
                samples(9,:) double
            end


            % determine particle charge
            switch options.type
                case 'p'
                    charge = 1;
                    mass = Consts.mp*Consts.c^2/(1e9*Consts.e);
                    fm = Consts.initialFourMomentum;
                    name = "proton";
                case 'e'
                    charge = -1;
                    mass = Consts.me*Consts.c^2/(1e9*Consts.e);
                    fm = [7000;0;0;sqrt(7000^2-mass^2)];
                    name = "electron";
                otherwise
                    charge = 1;
                    mass = Consts.mp*Consts.c^2/(1e9*Consts.e);
                    fm = Consts.initialFourMomentum;
                    name = "proton";
            end
            if n>1
                name = name + "s";
            end
            options.echo = ismember(options.echo,"on");
            if options.echo
                fprintf("Generating %u %s...\n",n,name)
            end
            % begin timer
            tic
            if isempty(obj.pdf)
                obj.pdf = @(~) 1;
            end
            cdfxdir = obj.cdfx;
            cdfydir = obj.cdfy;

            % perform checks on cdf table:
            % check table is sorted and monotonic(ish) increasing
            assert(all(all(cdfxdir(:,1:end-1)<=cdfxdir(:,2:end))))
            assert(all(all(cdfydir(:,1:end-1)<=cdfydir(:,2:end))))
            % check that cumulative probability is 1 at final value
            assert(cdfxdir(2,end)==1)
            assert(cdfydir(2,end)==1)

            % random sampling of x and y positions using inverse cdf table
            r = rand(2,n);
            xpos = interp1(cdfxdir(2,:),cdfxdir(1,:),r(1,:));
            ypos = interp1(cdfydir(2,:),cdfydir(1,:),r(2,:));
            % z positions at upstream edge 
            zpos = obj.zmin(ypos);

            
            samples = zeros(9,n);
            samples(1:3,:) = [xpos;ypos;zpos];
            samples(4:9,:) = repmat([fm;charge;mass],1,n);

%             fn = @(cdfxdir,cdfydir,dz,type) obj.generateSample(cdfxdir,cdfydir);
%             fn([0 0.5 1],[0 0.5 1],1,'p');
%             fn([0 0.5 1],[0 0.5 1],1,'p')
%             parfor i=1:n
%                 positions(:,i) = fn(cdfxdir,cdfydir,dz)
%             end
%             samples(1:3,:) = positions;
%             samples(4:9,:) = repmat([fm;charge;mass],1,n);
            tElapsed = toc;
            seconds = rem(tElapsed,60);
            minutes = (tElapsed-seconds)./60;
            % display simulation time and info
            if options.echo
                fprintf("Sample generation complete. Time taken was %d min %.5f s\n",minutes,seconds)
            end
        end

        function f = view(obj,f)
            if nargin==1
                f = figure;
            end
            geomViewStyle = "b";
            hold on
            plot3(obj.xv,obj.yv,obj.zvMin,geomViewStyle)
            plot3(obj.xv,obj.yv,obj.zvMax,geomViewStyle)
            view(3)
            % plot lines linking upper and lower polygons
            for i = 1:length(obj.xv)
                currLine = [obj.xv(i) obj.xv(i);obj.yv(i) obj.yv(i);obj.zvMin(i) obj.zvMax(i)];
                plot3(currLine(1,:),currLine(2,:),currLine(3,:),geomViewStyle);
            end
            xlabel("x")
            ylabel("y")
            zlabel("Trajectory direction - z")
            hold off
        end
    end

    methods % public, geometrical
        function b = isInside(obj,position)
            arguments(Input)
                obj(1,1) Geometry
                position(3,:) double    % [x;y;z] vals (m) to query
            end
            arguments(Output)
                b(1,:) logical          % logical for each position being inside the geometry
            end
            % check for separate volumes
            nanidx = [find(isnan(obj.xv)),length(obj.xv)+1];
            % split xv and yv at NaN values
            jlow = 1;
            b = zeros(1,size(position,2));
            for i = 1:length(nanidx)
                jhigh = nanidx(i)-1;
                b = b | isInsideSingle(position,obj.zvMin(jlow:jhigh),obj.zvMax(jlow:jhigh),obj.xv(jlow:jhigh),obj.yv(jlow:jhigh),obj.precision);
                jlow = nanidx(i)+1;
            end
        end

        function nextpos = edgeIntersect(obj,pos,dir)
            arguments(Input)
                obj(1,1) Geometry
                pos(3,:) double {mustBeFinite,Geometry.mustBeInside(pos,obj)}
                dir(3,:) double {mustBeFinite,Geometry.mustBeDirection(dir)}
            end
            arguments(Output)
                nextpos(3,:) double {mustBeFinite}
            end

            % check for separate volumes
            nanidx = [find(isnan(obj.xv)) length(obj.xv)+1];
            nextpos = zeros(3,size(pos,2));
            jlow = 1;
            for i=1:length(nanidx)
                jhigh = nanidx(i)-1;
                xvi = obj.xv(jlow:jhigh);
                yvi = obj.yv(jlow:jhigh);
                zvMaxi = obj.zvMax(jlow:jhigh);
                zvMini = obj.zvMin(jlow:jhigh);
                b = isInsideSingle(pos,zvMini,zvMaxi,xvi,yvi,obj.precision);
                %                 newpositions = obj.edgeIntersectSingle(pos(:,b),dir(:,b),xvi,yvi,zvMini,zvMaxi);
                nextpos(:,b) = obj.edgeIntersectSingle(pos(:,b),dir(:,b),xvi,yvi,zvMini,zvMaxi);
                jlow = nanidx(i)+1;
            end
        end
    end

    methods % geometry getters/setters
        function xmax = get.xmax(obj)
            xmax = max(obj.xv);
        end

        function xmin = get.xmin(obj)
            xmin = min(obj.xv);
        end

        function ymax = get.ymax(obj)
            ymax = max(obj.yv);
        end

        function ymin = get.ymin(obj)
            ymin = min(obj.yv);
        end

        function zmax = zmax(obj,y)
            % change so it works with xrotation
            arguments
                obj(1,1) Geometry
                y(1,:) double
            end
            yvals = [obj.ymin,obj.ymax];
            zvals = [min(obj.zvMax),max(obj.zvMax)];
            zmax = interp1(yvals,zvals,y);
        end

        function zmin = zmin(obj,y)
            arguments
                obj(1,1) Geometry
                y(1,:) double
            end
            yvals = [obj.ymin,obj.ymax];
            zvals = [min(obj.zvMin),max(obj.zvMin)];
            zmin = interp1(yvals,zvals,y);
        end

        function set.xRotation(obj,rotation)
            arguments
                obj(1,1) Geometry
                rotation(1,1) double {mustBeNonnegative,mustBeLessThan(rotation,90)}
            end
            obj.xRotation = rotation;
            obj.onAngleupdate();
        end
    
        function xlim = get.xLimits(obj)
            nanidx = [find(isnan(obj.xv)) length(obj.xv)+1];
            xlim = zeros(2,length(nanidx));
            jlow = 1;
            for i=1:length(nanidx)
                jhigh = nanidx(i)-1;
                maxX = max(obj.xv(jlow:jhigh));
                minX = min(obj.xv(jlow:jhigh));
                xlim(:,i) = [minX;maxX];
            end
        end

        function ylim = get.yLimits(obj)
            nanidx = [find(isnan(obj.yv)) length(obj.yv)+1];
            ylim = zeros(2,length(nanidx));
            jlow = 1;
            for i=1:length(nanidx)
                jhigh = nanidx(i)-1;
                maxY = max(obj.yv(jlow:jhigh));
                minY = min(obj.yv(jlow:jhigh));
                ylim(:,i) = [minY;maxY];
            end
        end
    end

    methods % simulation getters/setters
        function set.pdf(obj,newpdf)
            arguments
                obj(1,1) Geometry
                newpdf(1,1) function_handle
            end
            % validate pdf
            obj.pdf = newpdf;
        end

        function cdf = get.cdfx(obj)
            if isempty(obj.pdf)
                PDF = @(x) 1;
            else
                PDF = obj.pdf;
            end
            cdf = cdfgen(obj.xv,PDF);
        end

        function cdf = get.cdfy(obj)
            if isempty(obj.pdf)
                PDF = @(x) 1;
            else
                PDF = obj.pdf;
            end
            cdf = cdfgen(obj.yv,PDF);
        end
    end

    methods % to be overloaded if needed
        function [b,reEntryPos] = willReEnter(obj,pos,direction)
            arguments(Input)
                obj(1,1) Geometry
                pos(3,:) double
                direction(3,:) double
            end
            arguments(Output)
                b(1,:) logical
                reEntryPos(3,:) double
            end
            b = false(1,size(pos,2));
            reEntryPos = [];
        end
    end

    methods(Access=protected,Static) % input validation methods
        function mustBeInside(pos,obj)
            i = ~obj.isInside(pos);
            if any(i)
                eidType = "mustBeInside:notInside";
                msgType = "Not all samples are inside the geometry volume";
                throwAsCaller(MException(eidType,msgType))
            end
        end
        function mustBeSameLength(inp1,inp2)
            l1 = length(inp1);
            l2 = length(inp2);
            i = l1==l2;
            if ~i
                eidType = "mustBeSameLength:notSameLength";
                msgType = "Position and direction inputs are not the same size";
                throwAsCaller(MException(eidType,msgType))
            end
        end
        function mustBeDirection(inp)
            b = (inp(1,:)==0 & inp(2,:)==0 & inp(3,:)==0);
            if any(b)
                eidType = "mustBeDirection:zeroVector";
                msgType = "Some directions are the zero vector";
                throwAsCaller(MException(eidType,msgType))
            end
        end
    end

    methods(Access=private)

        function onAngleupdate(obj)
            % new angle in radians
            rad = obj.xRotation/180*pi;
            % change zv values
            zv = obj.yv*tan(rad);
            obj.zvMin = zv - min(zv); % ensure that at least one vertex is at z = 0
            obj.zvMax = obj.zvMin + obj.thickness;
            % check the geometrySize value again
            geometrySize = max(abs([obj.xv obj.yv obj.zvMin obj.zvMax]));
            obj.precision = abs(ceil(log10(eps(geometrySize))));
        end

        function nextpos = edgeIntersectSingle(obj,pos,dir,xv,yv,zvMin,zvMax)
            arguments(Output)
                nextpos(3,:) double {mustBeFinite}
            end
%             theta = obj.xRotation/180*pi;
%             edgeNormals = [0 0;sin(theta) -sin(theta);-cos(theta) cos(theta)];
%             edgeOriginOffsets = [dot([xv(1);yv(1);zvMin(1)],edgeNormals(:,1)),dot([xv(1);yv(1);zvMax(1)],edgeNormals(:,2))];
            aEdge = [xv(1) xv(1);yv(1) yv(1);zvMin(1) zvMax(1)];
            bEdge = [xv(2) xv(2);yv(2) yv(2);zvMin(2) zvMax(2)];
            cEdge = [xv(3) xv(3);yv(3) yv(3);zvMin(3) zvMax(3)];
            edgeNormals = cross(aEdge-bEdge,cEdge-bEdge);
            edgeNormals(:,1) = -edgeNormals(:,1);
            edgeNormals = edgeNormals./vecnorm(edgeNormals);
            edgeOriginOffsets = dot(edgeNormals,aEdge);
            assert(all(edgeOriginOffsets == dot(edgeNormals,bEdge)))
            
            % use 3 points on each face to find plane equation
            a = [xv(1:end-1);yv(1:end-1);zvMin(1:end-1)];
            b = [xv(1:end-1);yv(1:end-1);zvMax(1:end-1)];
            c = [xv(2:end);yv(2:end);zvMax(2:end)];
            % find normals for each plane
            polyNormals = cross(a-b,c-b);
            polyNormals = polyNormals./vecnorm(polyNormals);
            % find orogin offset for each plane
            polyOriginOffsets = dot(polyNormals,a,1);
            normals = [polyNormals edgeNormals];
            originOffsets = [polyOriginOffsets edgeOriginOffsets];
            % solve for t - nFaces x nPositions matrix
            pdotn = normals'*pos;
            vdotn = normals'*dir;
            t = -(pdotn-repmat(originOffsets',1,size(pos,2)))./vdotn;
            % fix small numerical errors
            t = round(t,obj.precision);
            % find the least positive value of t in each column
            chosenTVals = minPositive(t);

            % plug into expression for projected position
            % again fix numerical errors - this time to eps of position
            posPrecision = min(abs(ceil(log10(eps(pos)))))-1;
            dPos = round(dir.*chosenTVals,min(obj.precision,posPrecision));
            nextpos = pos + dPos;

            % z position numerical error bodge
%             if nextpos(3) > obj.zmax(nextpos(2))
%                 nextpos(3) = obj.zmax(nextpos(2));
%             end
            
%             if any(~obj.isInside(nextpos))
%                 throwAsCaller(MException("edgeIntersectSingle:projectedOutside","Projected edge intersection is outside the geometry"))
%             end
%             Geometry.mustBeInside(nextpos,obj)
        end
    end
end

function b = isInsideSingle(pos,zvMin,zvMax,xv,yv,precision)
x = pos(1,:);
y = pos(2,:);
z = round(pos(3,:),precision);
yBounds = [min(yv),max(yv)];
% if ~isreal(y)
%     y
% end
zMax = round(interp1(yBounds,[min(zvMax),max(zvMax)],y),precision);
zMin = round(interp1(yBounds,[min(zvMin),max(zvMin)],y),precision);
b = inpolygon(x,y,xv,yv) & (z<=zMax & z>=zMin);
end

function v = minPositive(mat)
mat(mat<=0) = Inf;
% output vector allocation
v = zeros(1,size(mat,2));
% columns where there is at least one strictly positive value
positiveIdx = ~all(mat==Inf);
% find min only of columns with at leat one strictly positive value
v(positiveIdx) = min(mat(:,positiveIdx));
end

function [cdf,totalProbability] = cdfgen(vertices,pdf)
nanidx = [find(isnan(vertices)) length(vertices) + 1];
nSections = length(nanidx);
sectionLength = nanidx(1)-1;
sections = zeros(nSections,sectionLength);
jlow = 1;
for i = 1:nSections
    jhigh = nanidx(i)-1;
    sections(i,:) = vertices(jlow:jhigh);
    jlow = nanidx(i) + 1;
end
% we now have a matrix with rows that are separate x sections
sections = unique(sections,'rows');
% recalculate the number of sections
nSections = size(sections,1);
% get resolution from runtimedata
resolution = Consts.data.cdfResolution;
% preallocate output cdf array
cdf = zeros(2,resolution*nSections); %
% preallocate xvals for integration of each section, and cdfvals
% resulting from integration along the xvals of each section
cdfvals = zeros(nSections,resolution);
xvals = zeros(nSections,resolution); % rows are xvals for section i
fun = Consts.data.gridDensityFunction;
factor = Consts.data.gridDensityMultiplier;
for i = 1:nSections
    % calculate the x values for each section
    xmax = max(sections(i,:));
    xmin = min(sections(i,:));
    uniformgrid = linspace(xmin,xmax,resolution);
    xvals(i,:) = nonUniformGrid(uniformgrid,factor,fun);
    for j = 2:length(xvals(i,:))
        % integrate between j-1 and jth elements of xval
        dx = xvals(i,j)-xvals(i,j-1);
        p = 0.5*(pdf(xvals(i,j-1))+pdf(xvals(i,j)));
        dP = dx*p;
        cdfvals(i,j) = cdfvals(i,j-1) + dP;
    end
end
% reorder xvals (and cdfvals) so that x value increases down the rows
[~,I] = sort(max(xvals,[],2));
xvals = xvals(I,:);
cdfvals = cdfvals(I,:);
% make each row increase by the last entry of the previous rows (cumulative)
if nSections>1
    for i = 2:nSections
        cdfvals(i,:) = cdfvals(i,:) + cdfvals(i-1,end);
    end
end
cdf(1,:) = reshape(xvals',[],nSections*resolution);
cdf(2,:) = reshape(cdfvals',[],nSections*resolution);
totalProbability = cdf(2,end);
if totalProbability==0
    eid = "cdfvals:zeroProbability";
    msg = "The pdf provided is too small over the geometry and integrates to zero";
    throw(MException(eid,msg));
end
% normalise
cdf(2,:) = cdf(2,:)./cdf(2,end);
% make sure cdf is monotonically increasing
Inew = (1:nSections-1)*resolution + 1;
cdf(2,Inew) = cdf(2,Inew) + eps(cdf(2,Inew));
% clean out repeating ones and zeros from the start and end of cdf
nMinus = floor(length(cdf(2,:))/2);
nPlus = nMinus+1;
[Pminus,iaMinus,~] = unique(cdf(2,1:nMinus),'first');
xMinus = cdf(1,1:nMinus);
[Pplus,iaPlus,~] = unique(cdf(2,nPlus:end),'last');
xPlus = cdf(1,nPlus:end);
cdf = [xMinus(iaMinus) xPlus(iaPlus);Pminus Pplus];
end
