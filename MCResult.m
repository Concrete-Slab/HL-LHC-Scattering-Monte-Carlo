classdef MCResult
    %MCRESULT Summary of this class goes here
    %   Detailed explanation goes here

    properties(SetAccess=private)
        primaries(9,:) double
        primaryHistory(1,:) ParticleHistory
        secondaries(9,:) double
        secondaryHistory(1,:) ParticleHistory
    end

    properties(SetAccess=private) % NaN separated particle tracks
        primaryTracks(3,:) double
        secondaryTracks(3,:) double
    end

    properties(SetAccess=private) % energy deposition variables
        projectedPrimaries(9,:) double
        projectedSecondaries(9,:) double
    end

    properties(Access=private)
        mcinput(1,:) MCInput = MCInput.empty(1,0);
    end

    properties(Dependent,SetAccess=private)
        input(9,:) double
        geometry(1,1) Geometry
        material(1,1) Material
        projectedParticles(9,:) double
    end

    methods % constructor
        function obj = MCResult(primaries,secondarySamples,initial,primaryHistory)
            %MCRESULT Construct an instance of this class
            %   Detailed explanation goes here
            arguments(Input)
                primaries(9,:) double {mustBeFinite} = double.empty(9,0)
                secondarySamples(1,:) SecondarySamples = SecondarySamples.empty
                initial(1,1) MCInput = MCInput
                primaryHistory(1,:) ParticleHistory = ParticleHistory.empty
            end

            % assign primary properties and input property
            obj.primaries = primaries;
            obj.mcinput = initial;
            obj.primaryHistory = primaryHistory;

            % ascertain how many secondaries were produced
            nSecondaries = 0;
            for i = 1:length(secondarySamples)
                nSecondaries = nSecondaries + secondarySamples(i).n;
            end

            % now preallocate the secondary sample and history arrays
            secs = zeros(9,nSecondaries);
            hists = repmat(ParticleHistory,1,nSecondaries);

            % fill the arrays
            j = 1;
            for i = 1:length(secondarySamples)
                sizeToAdd = secondarySamples(i).n;
                jnext = j + sizeToAdd-1;
                secs(:,j:jnext) = secondarySamples(i).secondaries;
                hists(:,j:jnext) = secondarySamples(i).secondaryHistories;
                j = jnext+1;
            end
            obj.secondaries = secs;

            % generate particle history object corresponding to each secondary
            obj.secondaryHistory = hists;

            % generate NaN-separated secondary tracks
            obj.secondaryTracks = getTracks(obj.secondaryHistory);

            % generate NaN-separated primary tracks (history may not exist)
            if ~isempty(obj.primaryHistory)
                obj.primaryTracks = getTracks(obj.primaryHistory);
            else
                obj.primaryTracks = zeros(3,3*size(obj.primaries,2));
                obj.primaryTracks(:,1:3:end) = obj.input(1:3,:);
                obj.primaryTracks(:,2:3:end) = obj.primaries(1:3,:);
                obj.primaryTracks(:,3:3:end) = NaN;
            end

            % calculate the extrapolated position of impact with aperture
            % assumes linear motion

            R = Consts.TubeRadius;
            % formula in terms of transverse displacement only from origin,
            % quadratic for scalar in direction of momentum
            projectedParticles = [obj.primaries,obj.secondaries];
            % starting positions of all particles
            initialPos = projectedParticles(1:3,:);
            % directions don't have to be normalised, only relative values matter
            initialDirection = [obj.primaries(5:7,:),obj.secondaries(5:7,:)];
            % indices of particles that will eventually hit tube - not
            % parallel to z axis, and above minimum kinetic energy
            mke = Consts.MinimumKineticEnergy;
            primaryParticlesIndex = ~(obj.primaries(5,:)==0 & obj.primaries(6,:)==0) & (obj.primaries(4,:)-obj.primaries(9,:))>mke;
            secondaryParticlesIndex = [false(1,length(primaryParticlesIndex)), ~(obj.secondaries(5,:)==0 & obj.secondaries(6,:)==0) & (obj.secondaries(4,:)-obj.secondaries(9,:))>mke];
            % dot product of each position column with itself
            xdotx = dot(initialPos(1:2,:),initialPos(1:2,:));
            % dot product of each position and corresponding direction
            xdotp = dot(initialPos(1:2,:),initialDirection(1:2,:));
            pdotp = dot(initialDirection(1:2,:),initialDirection(1:2,:));
            % take positive root of each calculation
            alpha = (-xdotp+sqrt(xdotp.^2-(pdotp.*(xdotx-R^2))))./pdotp;
            alpha = repmat(alpha,3,1);
            % final position = x + alpha*p
            projectedParticles(1:3,:) = initialPos + alpha.*initialDirection;
            % remove particles that never reach z axis
            obj.projectedPrimaries = projectedParticles(:,primaryParticlesIndex);
            obj.projectedSecondaries = projectedParticles(:,secondaryParticlesIndex);
%             obj.projectedParticles = obj.projectedParticles(:,projectedParticlesIndex);
            
%             obj.projectedParticles = obj.projectedParticles(:,isreal(obj.projectedParticles));
        end
    end

    methods % getters
        function g = get.geometry(obj)
            if isempty(obj.mcinput)
                g = SymRectangles.empty(1,0);
            else
                g = obj.mcinput.geometry;
            end
        end
        function m = get.material(obj)
            if isempty(obj.mcinput)
                m = Material.empty(1,0);
            else
                m = obj.mcinput.material;
            end
        end
        function s = get.input(obj)
            if isempty(obj.mcinput)
                s = double.empty(8,0);
            else
                s = obj.mcinput.samples;
            end
        end
        function pp = get.projectedParticles(obj)
            pp = [obj.projectedPrimaries,obj.projectedSecondaries];
        end
%         function pp = get.projectedPrimaries(obj)
%             I = find(~(obj.primaries(5,:)==0 & obj.primaries(6,:)==0));
%             pp = obj.projectedParticles(:,1:size(obj.primaries(:,I),2));
%         end
%         function ps = get.projectedSecondaries(obj)
%             ps = obj.projectedParticles(:,(size(obj.primaries,2)+1):end);
%         end
    end

    methods % public, visualisation
        function viewTracks(obj)
            %VIEWTRACKS Visualises the geometry along with the primary
            %particle tracks in purple and secondary tracks in green
            obj.mcinput.view
            hold on
            pt = obj.primaryTracks;
            st = obj.secondaryTracks;
            plot3(pt(1,:),pt(2,:),pt(3,:),'m')
            if ~isempty(obj.secondaries)
                plot3(st(1,:),st(2,:),st(3,:),'g')
            end
            hold off
        end

        function viewCompleteTracks(obj)
            %VIEWCOMPLETETRACKS Extrapolates the particle tracks to the
            %intersection with the beam tube, displays with geometry
            obj.mcinput.view;
            pt = obj.primaryTracks;
            st = obj.secondaryTracks;
            % replace the final position in each track with the
            % extrapolated position.
            % These column indices are immediately before NaN separators
            primFinalPos = find(isnan(pt(1,:)))-1;
            secFinalPos = find(isnan(st(1,:)))-1;
            pt(:,primFinalPos) = obj.projectedPrimaries(1:3,:);
            if isempty(obj.projectedSecondaries) % 0 secondaries
                st = double.empty(3,0);
            elseif isempty(secFinalPos)         % 1 secondary
                st = [obj.secondaryTracks(1:3,1),obj.projectedSecondaries(1:3,1)];
            else                                % more than 1 secondary
                st(:,secFinalPos) = obj.projectedSecondaries(1:3,:);
            end
            hold on
            plot3(pt(1,:),pt(2,:),pt(3,:),'m')
            if ~isempty(obj.secondaries)
                plot3(st(1,:),st(2,:),st(3,:),'g')
            end
            theta = linspace(-pi,pi,100);
            x = Consts.TubeRadius*sin(theta);
            y = Consts.TubeRadius*cos(theta);
            z = obj.mcinput.geometry.thickness*ones(1,length(theta));
            plot3(x,y,z)
            set(gca,'ZScale','log')
            hold off
        end

        function viewChanges(obj)
            %VIEWCHANGES Plots the change in key variables between the
            %particle entering and leaving the screen
            dSample = obj.mcinput.samples-obj.primaries;
            viewChangesDisplayStyle = 'r.';
            figure
            tiledlayout(3,3)
            nexttile
            plot(dSample(1,:),dSample(2,:),viewChangesDisplayStyle)
            xlabel("\Delta{x}, m")
            ylabel("\Delta{y}, m")
            nexttile
            plot(dSample(5,:),dSample(6,:),viewChangesDisplayStyle)
            xlabel("\Delta{p_x}, GeV/c")
            ylabel("\Delta{p_y}, GeV/c")
            nexttile
            plot(dSample(1,:),dSample(5,:),viewChangesDisplayStyle)
            xlabel("\Delta{x}, m")
            ylabel("\Delta{p_x}, GeV/c")
            nexttile
            plot(dSample(2,:),dSample(6,:),viewChangesDisplayStyle)
            xlabel("\Delta{y}, m")
            ylabel("\Delta{p_y}, GeV/c")
        end

        function viewIntersections(obj,options)
            arguments
                obj(1,1) MCResult
                options.scalefactor(1,1) double {mustBePositive} = 1;
                options.MagnetMaterial(1,1) Material = Consts.MagnetMaterial
                options.zGrid(1,:) double = []% m
                options.yGrid(1,:) double = []% m
                options.MaxZPoints(1,:) double {mustBeInteger,mustBeScalarOrEmpty} = []
                options.MaxThetaPoints(1,:) double {mustBeInteger,mustBeScalarOrEmpty} = []
                options.echo(1,1) string {mustBeMember(options.echo,["on","off"])} = "on";
            end
            thetaPrimary = Consts.TubeRadius.*atan2(obj.projectedPrimaries(2,:),obj.projectedPrimaries(1,:));
            zPrimary = obj.projectedPrimaries(3,:);
            thetaSecondary = Consts.TubeRadius.*atan2(obj.projectedSecondaries(2,:),obj.projectedSecondaries(1,:));
            zSecondary = obj.projectedSecondaries(3,:);
            tiledlayout(2,1)
            nexttile;
            % make axes on both tiles equal somehow...
            semilogx(zPrimary,thetaPrimary,"m.")
            hold on
            semilogx(zSecondary,thetaSecondary,"g.")
            xlabel("Distance along beam tube, z")
            ylabel("Circumferential distance, R\theta")
            hold off
            ax = nexttile;
            [Z,Y,E] = obj.energyDistribution("Material",options.MagnetMaterial,"zGrid",options.zGrid,"yGrid",options.yGrid,"MaxZPoints",options.MaxZPoints,"MaxThetaPoints",options.MaxThetaPoints,"echo",options.echo);
            % E is in J/m^3 - change to W/m^3 and then to mw/cm^3
            E = E.*options.scalefactor.*1e-3;
            surf(Z,Y,E)
            set(ax,'xscale','log')
            xlabel("Distance along beam tube, m")
            ylabel("Distance around circumference, m")
            zlabel("Volumetric energy deposition, mW {cm}^{-3}")
        end
    end

    methods % public, data
        function [Z,Y,E] = energyDistribution(obj,options)
            arguments(Input)
                obj(1,1) MCResult
                options.Material(1,1) Material = Consts.MagnetMaterial;
                options.zGrid(1,:) double =[] % m
                options.yGrid(1,:) double =[] % m
                options.MaxZPoints(1,:) double {mustBeInteger,mustBeScalarOrEmpty} = []
                options.MaxThetaPoints(1,:) double {mustBeInteger,mustBeScalarOrEmpty} = []
                options.echo(1,1) string {mustBeMember(options.echo,["on","off"])} = "on";
            end
            arguments(Output)
                Z(:,:) double % grid, m
                Y(:,:) double % grid, m
                E(:,:) double % grid, GeV/m^3
            end
            options.echo = ismember(options.echo,"on");
            % input parsing
            if isempty(options.MaxThetaPoints)
                options.MaxThetaPoints = options.MaxZPoints;
            end

            % energy deposition of each particle
            gamma = obj.projectedParticles(4,:)./obj.projectedParticles(9,:);
            beta = sqrt(1-1./gamma.^2);
            % energy deposition in GeV/m
            dEdx = options.Material.dEdx(beta,obj.projectedParticles(9,:),obj.projectedParticles(8,:));
            
            % total energy deposition
            totE = sum(dEdx);

            % positions of samples on z-theta*R graph
            z = obj.projectedParticles(3,:);
            theta = sign(obj.projectedParticles(2,:)).*acos(obj.projectedParticles(1,:)./sqrt(obj.projectedParticles(1,:).^2+obj.projectedParticles(2,:).^2));
            y = theta*Consts.TubeRadius;
            % store zysamples on the GPU to speed it up for larger samples
            zySamples = gpuArray([z;y])';

            % grid initialisation
            if isempty(options.zGrid)
                zGrid = equiprobableGrid(z,options.MaxZPoints);
            else
                zGrid = options.zGrid;
            end
            if isempty(options.yGrid)
                yGrid = equiprobableGrid(y,options.MaxThetaPoints);
            else
                yGrid = options.yGrid;
            end
            [Z,Y] = meshgrid(zGrid,yGrid);
            if options.echo
                % begin performance timer
                fprintf("Performing kernel energy density estimation onto %d x %d grid\n",size(Z,1),size(Z,2));
                tic
            end
            % KDE grid input (must be nx2), stored on GPU
            zyGrid = gpuArray([Z(:) Y(:)]);
            % KDE sample input is zySamples from before

            % define weights for the kernel density estimation, GPU array
            energyWeights = gpuArray(dEdx./totE);
            % perform a weighted KDE calculation
            [weightedPdf,~] = ksdensity(zySamples,zyGrid,'Weights',energyWeights,'Function','pdf');
            % multiply pdf by total energy deposition to get energy
            % deposition density stored on GPU
            Egpu = reshape(weightedPdf.*totE,size(Z));
            % energy deposition (main memory) in GeV/m^3
            E = gather(Egpu);
            % energy deposition in J/m^3
            E = E*1e9*Consts.e;
            if options.echo
                tElapsed = toc;
                seconds = rem(tElapsed,60);
                minutes = (tElapsed-seconds)./60;
                fprintf("Energy deposition calculated in %d min %.5f\n",minutes,seconds)
            end
        end

        function [Z,Y,SV] = protonDensity(obj,protonFlux,options)
            arguments(Input)
                obj(1,1) MCResult
                protonFlux(1,1) double % number of protons per second on screen
                options.Material(1,1) Material = "Cu";
                options.MaxZPoints(1,:) double {mustBeInteger,mustBeScalarOrEmpty} = []
                options.MaxThetaPoints(1,:) double {mustBeInteger,mustBeScalarOrEmpty} = []
                options.echo(1,1) string {mustBeMember(options.echo,["on","off"])} = "on";
            end
            arguments(Output)
                Z(:,:) double % grid, m
                Y(:,:) double % grid, m
                SV(:,:) double % grid, GeV/m^3
            end
            options.echo = ismember(options.echo,"on");
            % input parsing
            if isempty(options.MaxThetaPoints)
                options.MaxThetaPoints = options.MaxZPoints;
            end

            % positions of samples on z-theta*R graph
            zPrim = obj.projectedPrimaries(3,:);
            zSec = obj.projectedSecondaries(3,:);
            thetaPrim = sign(obj.projectedPrimaries(2,:)).*acos(obj.projectedPrimaries(1,:)./sqrt(obj.projectedPrimaries(1,:).^2+obj.projectedPrimaries(2,:).^2));
            thetaSec = sign(obj.projectedSecondaries(2,:)).*acos(obj.projectedSecondaries(1,:)./sqrt(obj.projectedSecondaries(1,:).^2+obj.projectedSecondaries(2,:).^2));
            yPrim = thetaPrim*Consts.TubeRadius;
            ySec = thetaSec*Consts.TubeRadius;
            % store zysamples on the GPU to speed it up for larger samples
            zySamplesPrim = gpuArray([zPrim;yPrim])';
            zySamplesSec = gpuArray([zSec;ySec])';
            % grid initialisation
            zGridPrim = equiprobableGrid(zPrim,options.MaxZPoints);
            zGridSec = equiprobableGrid(zSec,options.MaxZPoints);
            % keep yGrid the same for both primaries and secondaries
            yGrid = equiprobableGrid(y,options.MaxThetaPoints);
            
            [ZPrim,Y] = meshgrid(zGridPrim,yGrid);
            [ZSec,~] = meshgrid(zGridSec,yGrid);
            if options.echo
                % begin performance timer
                fprintf("Performing kernel energy density estimation onto %d x %d grid\n",size(Z,1),size(Z,2));
                tic
            end
            % KDE grid input (must be nx2), stored on GPU
            zyGridPrim = gpuArray([ZPrim(:) Y(:)]);
            zyGridSec = gpuArray([ZSec(:), Y(:)]);
            % KDE sample input is zySamples from before
            [PDFPrim,~] = ksdensity(zySamplesPrim,zyGridPrim,'Function','pdf');
            if ~isempty(zySamplesSec)
                [PDFSec,~] = ksdensity(zySamplesSec,zyGridSec);
            else
                PDFSec = double.empty(size(Y,2),0);
            end
            propPrimariesIntersecting = size(obj.projectedPrimaries,2)./size(obj.primaries,2);
            Z = [ZSec ZPrim];
            PDF = [PDFSec PDFPrim];
            SV = PDF.*propPrimariesIntersecting.*protonFlux;
        end
    end
end

function tr = getTracks(history)
nTracks = length(history);
nPositions = 0;
for i = 1:nTracks
    nPositions = nPositions + size(history(i).positions,2);
end
% preallocate - add nTracks to account for each NaN column
tr = zeros(3,nPositions + nTracks);
j = 1;
for i = 1:nTracks
    lengthToAdd = history(i).n;
    jnext = j + lengthToAdd - 1;
    tr(:,j:jnext) = history(i).positions;
    tr(:,jnext+1) = NaN;
    j = jnext + 2;
end
end

function [gridVals,gridSteps] = equiprobableGrid(vec,maxK)
arguments(Input)
    vec(1,:) double % all values of a certain variable
    maxK(1,:) double {mustBeScalarOrEmpty,mustBeInteger,mustBePositive} = []; % maximum number of grid points
end
% vector needs to be sorted before starting
vec = sort(vec);
n = length(vec);
if n==0
    gridVals = [];
    gridSteps = [];
    return
end
if isempty(maxK)
    k = ceil(min(2*n^(2/5),n/5)); % NIST formula, at least 5 samples per division
else
    k = min(ceil(min(2*n^(2/5),n/5)),maxK);
end
valsPerBin = round(n/k);
I = [1:valsPerBin:n n];
% mean grid positions
gridEdges = vec(I);
gridVals = 0.5*(gridEdges(1:end-1)+gridEdges(2:end));
gridSteps = gridEdges(2:end)-gridEdges(1:end-1);
end