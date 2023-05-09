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
            mkePrimary = Consts.MinimumKineticEnergy(obj.primaries(9,:),obj.primaries(8,:));
            mkeSecondary = Consts.MinimumKineticEnergy(obj.secondaries(9,:),obj.secondaries(8,:));
            primaryParticlesIndex = ~(obj.primaries(5,:)==0 & obj.primaries(6,:)==0) & (obj.primaries(4,:)-obj.primaries(9,:))>mkePrimary;
            secondaryParticlesIndex = [false(1,length(primaryParticlesIndex)), ~(obj.secondaries(5,:)==0 & obj.secondaries(6,:)==0) & (obj.secondaries(4,:)-obj.secondaries(9,:))>mkeSecondary];
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
            tiledlayout(2,2)
            nexttile
            plot(dSample(1,:),dSample(2,:),viewChangesDisplayStyle)
            xlabel("\Delta{x}, m")
            ylabel("\Delta{y}, m")
            nexttile
            plot(dSample(1,:),dSample(5,:),viewChangesDisplayStyle)
            xlabel("\Delta{x}, m")
            ylabel("\Delta{p_x}, GeV/c")
            nexttile
            plot(dSample(5,:),dSample(6,:),viewChangesDisplayStyle)
            xlabel("\Delta{p_x}, GeV/c")
            ylabel("\Delta{p_y}, GeV/c")
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
                options.kmeans(1,1) double {mustBeInteger} = 3;
                options.ShieldSecondaries(1,1) logical = 0;
                options.ShieldingMaterial(1,1) Material = "W"
                options.ShieldingDepth(1,1) double {mustBePositive} = Consts.TungstenDepth;
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
            [Z,Y,E] = obj.energyDistribution("Material",options.MagnetMaterial,"zGrid",options.zGrid,"yGrid",options.yGrid,"MaxZPoints",options.MaxZPoints,"MaxThetaPoints",options.MaxThetaPoints,"echo",options.echo,"kmeans",options.kmeans,"ShieldSecondaries",options.ShieldSecondaries,"ShieldingMaterial",options.ShieldingMaterial,"ShieldingDepth",options.ShieldingDepth);
            % E is in J/m^3 - change to W/m^3 and then to mw/cm^3
            E = E.*options.scalefactor.*1e-3;
            surf(Z,Y,E)
            set(ax,'xscale','log')
            xlabel("Distance along beam tube, m")
            ylabel("Distance around circumference, m")
            zlabel("Volumetric energy deposition, mW {cm}^{-3}")
        end

        function viewCSDARange(obj,material,options)
            arguments(Input)
                obj(1,1) MCResult
                material(1,1) Material = "W"
                options.dims(1,1) double {mustBeMember(options.dims,[2,3])} % number of plot axes
                options.include(1,1) string {mustBeMember(options.include,["primaries","secondaries","all"])} = "all"
                options.method(1,1) string {mustBeMember(options.method,["CSDA","True","Perpendicular"])} = "CSDA";
            end
            if options.method == "CSDA"
                [delXPrimaries,delXSecondaries] = obj.getCSDARange("include",options.include);
            else
                [delXPrimaries,delXSecondaries] = obj.getRange("include",options.include);
                if options.method=="Perpendicular"
                    if options.include == "primaries"||options.include=="all"
                        primaryDirections = sym(obj.projectedPrimaries(5:7,:));
                        primarySines = double(sin(acos(primaryDirections(3,:)./sqrt(primaryDirections(3,:).^2+primaryDirections(2,:).^2+primaryDirections(1,:).^2))));
                        delXPrimaries = delXPrimaries.*primarySines;
                    end
                    if options.include=="secondaries"||options.include=="all"
                        secondaryDirections = sym(obj.projectedSecondaries(5:7,:));
                        secondarySines = double(sin(acos(secondaryDirections(3,:)./sqrt(secondaryDirections(3,:).^2+secondaryDirections(2,:).^2+secondaryDirections(1,:).^2))));
                        delXSecondaries = delXSecondaries.*secondarySines;
                    end
                end
            end
            
            if isempty(delXPrimaries)&&isempty(delXSecondaries)
                warning("No particles included")
                return
            end
            % plot results
            figure
            tiledlayout(1,2)
            nexttile
            if ~isempty(delXPrimaries)
                loglog(obj.projectedPrimaries(3,:),delXPrimaries,"m.");
                hold on
            end
            loglog(obj.projectedSecondaries(3,:),delXSecondaries,"g.");
            hold off
            xlabel("Distance along beam axis z, m")
            strLabel = sprintf("%s range in %s, m",options.method,material.name);
            ylabel(strLabel);
            if options.include == "all"
                legend("Secondaries","Primaries");
            end
            nexttile
            if ~isempty(delXPrimaries)
                loglog(obj.projectedPrimaries(3,:),obj.projectedPrimaries(4,:),"m.");
                hold on
            end
            loglog(obj.projectedSecondaries(3,:),obj.projectedSecondaries(4,:),"g.");
            hold off
            xlabel("Distance along beam axis z, m")
            ylabel("Final particle energy, GeV")
            if options.include =="all"
                legend("Secondaries","Primaries");
            end
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
                options.ShieldSecondaries(1,1) logical = 0;
                options.ShieldingMaterial(1,1) Material = "W"
                options.ShieldingDepth(1,1) double {mustBePositive} = Consts.TungstenDepth
                options.echo(1,1) string {mustBeMember(options.echo,["on","off"])} = "on";
                options.kmeans(1,1) double {mustBeInteger} = 3
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
            if options.ShieldSecondaries
                if options.echo
                    fprintf("Calculating secondary ranges...")
                end
                % need to first filter out the secondaries without range
                [~,secondaryRanges] = obj.getRange(options.ShieldingMaterial,"include","secondaries","method","Perpendicular");
                includeIdx = secondaryRanges>options.ShieldingDepth;
                includedSecondaries = obj.projectedSecondaries(:,includeIdx);
%                 includedSecondaries = obj.projectedSecondaries;
                % now decrease their energy according to path length
                shieldingMaterial = options.ShieldingMaterial;
                includedSecondaries(4,:) = shieldingMaterial.energyChange(includedSecondaries(4,:)-includedSecondaries(9,:),includedSecondaries(9,:),includedSecondaries(8,:),secondaryRanges(includeIdx));
                % alter the momentum magnitude according to new energy
                includedSecondaries(5:7,:) = includedSecondaries(5:7,:).* sqrt((includedSecondaries(4,:).^2-includedSecondaries(9,:).^2)./(includedSecondaries(5,:).^2+includedSecondaries(6,:).^2+includedSecondaries(7,:).^2));
                % save secondaries
                projectedData = [obj.projectedPrimaries includedSecondaries];
            else
                projectedData = obj.projectedParticles;
            end

            beta = sqrt(1-1./(projectedData(4,:)./projectedData(9,:)).^2);
            % energy deposition in GeV/m
            dEdx = options.Material.dEdx(beta,projectedData(9,:),projectedData(8,:));
            % positions of samples on z-theta*R graph
            z = projectedData(3,:);
            theta = sign(projectedData(2,:)).*acos(projectedData(1,:)./sqrt(projectedData(1,:).^2+projectedData(2,:).^2));
            y = theta*Consts.TubeRadius;
            % store zysamples on the GPU to speed it up for larger samples
            zySamples = [z;y]';

            if isempty(options.yGrid)
                yGrid = equiprobableGrid(y,options.MaxThetaPoints);
            else
                yGrid = options.yGrid;
            end
            if options.echo
                % begin performance timer
                fprintf("Performing kernel energy density estimation");
                tic
            end
            % perform 3-mean clustering on log vals of z axis
            logZ = log10(abs(z));
            % TODO make k and optional parameter
            k = options.kmeans;
            if k==1
                clusterIdx = ones(1,length(z));
            else
                clusterIdx = kmeans(logZ',k)';
            end

            % not necessary for results, but helps for visualisation:
            % sort the indices so that the smallest mean is group 1 etc.
            firstZpositions = zeros(1,k);
            for i = 1:k
                zK = z(clusterIdx==i);
                firstZpositions(i) = zK(1);
            end
            % get sorting order
            [~,clusterOrder] = sort(firstZpositions);
            % sort the indices
            unsortedIdx = clusterIdx;
            for i = 1:k
                clusterIdx(unsortedIdx==clusterOrder(i)) = i;
            end

            % make a grid for each cluster, splice together for final grid
            nZGrid = 0;
            nj = zeros(1,k);
            for j = 1:k
                nj(j) = length(z(clusterIdx==j));
                if isempty(options.MaxZPoints)
                    nZGrid = nZGrid + ceil(min(2*nj(j)^(2/5),nj(j)/5));
                else
                    nZGrid = nZGrid + min(ceil(min(2*nj(j)^(2/5),nj(j)/5)),options.MaxZPoints);
                end
            end
            zGrid = zeros(1,nZGrid);

            cellSamples = cell(1,k);
            cellWeights = cell(1,k);
            currIdx = 1;
            for j = 1:k
                zK = z(clusterIdx==j);
                thisGrid = equiprobableGrid(zK,options.MaxZPoints);

                zGrid(currIdx:currIdx+length(thisGrid)-1) = thisGrid;
                currIdx = currIdx+length(thisGrid);

                cellSamples(j) = {zySamples(clusterIdx==j,:)};
                cellWeights(j) = {dEdx(clusterIdx==j)};
            end
            [Z,Y] = meshgrid(zGrid,yGrid);
            zyGrid = gpuArray([Z(:),Y(:)]);
            % perform kde on each cluster:
            Ek = gpuArray(zeros(size(zyGrid,1),k));
            parfor i = 1:k
                % find total energy of the cluster
                clusterdEdx = cellWeights{i};
                clusterEnergy = sum(clusterdEdx);
                % make weights from relative energy to cluster total
                clusterWeights = gpuArray(clusterdEdx./clusterEnergy)';
                % apply weighted kde to global grid
                clusterSamples = gpuArray(cellSamples{i});
                clusterPDF = ksdensity(clusterSamples,zyGrid,"weights",clusterWeights);
                % get cluster's contribution to total energy on grid
                Ek(:,i) = clusterEnergy*clusterPDF;
            end

            % now add the results of all 3 kde's together

            Evec = sum(Ek,2);

            % reshape the final result back to grid

            E = gather(reshape(Evec,size(Z)));








            %             % perform a weighted KDE calculation
            %             weightedPdf = ksdensity(zySamples,zyGrid,'Weights',energyWeights);
            %             % multiply pdf by total energy deposition to get energy
            %             % deposition density stored on GPU
            %             Egpu = reshape(weightedPdf.*totE,size(Z));
            %             % energy deposition (main memory) in GeV/m^3
            %             E = Egpu;
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
            yGrid = equiprobableGrid(yPrim,options.MaxThetaPoints);

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

        function [primaryCSDA,secondaryCSDA] = getCSDARange(obj,material,options)
            arguments(Input)
                obj(1,1) MCResult
                material(1,1) Material = "W"
                options.include(1,1) string {mustBeMember(options.include,["primaries","secondaries","all"])} = "all"
                
            end
            primaryCSDA = [];
            secondaryCSDA = [];
            f = @(ek,m,z) material.csdaRange(ek,m,z);
            if options.include == "primaries" ||options.include=="all"
                eKprimaries= obj.projectedPrimaries(4,:)-obj.projectedPrimaries(9,:);
                massPrimaries = obj.projectedPrimaries(9,:);
                chargePrimaries = obj.projectedPrimaries(8,:);
                primaryCSDA = zeros(1,length(massPrimaries));
                parfor i = 1:length(massPrimaries)
                    primaryCSDA(i) = f(eKprimaries(i),massPrimaries(i),chargePrimaries(i));
                end
            end

            if options.include == "secondaries"||options.include=="all"
                eKsecondaries = obj.projectedSecondaries(4,:)-obj.projectedSecondaries(9,:);
                massSecondaries = obj.projectedSecondaries(9,:);
                chargeSecondaries = obj.projectedSecondaries(8,:);
                secondaryCSDA = zeros(1,length(massSecondaries));
                parfor i = 1:length(massSecondaries)
                    secondaryCSDA(i) = f(eKsecondaries(i),massSecondaries(i),chargeSecondaries(i));
                end
            end
        end

        function E = maxEnergy(obj,options)
            arguments(Input)
                obj(1,1) MCResult
                options.Material(1,1) Material = Consts.MagnetMaterial;
                options.MaxZPoints(1,:) double {mustBeInteger,mustBeScalarOrEmpty} = []
                options.MaxThetaPoints(1,:) double {mustBeInteger,mustBeScalarOrEmpty} = []
                options.echo(1,1) string {mustBeMember(options.echo,["on","off"])} = "on";
                options.kmeans(1,1) double {mustBeInteger} = 3
            end
            arguments(Output)
                E(:,:) double % max, J/m^3
            end
            options.echo = ismember(options.echo,"on");
            % input parsing
            if isempty(options.MaxThetaPoints)
                options.MaxThetaPoints = options.MaxZPoints;
            end
            projectedData = obj.projectedParticles;
            beta = sqrt(1-1./(projectedData(4,:)./projectedData(9,:)).^2);
            % energy deposition in GeV/m
            dEdx = options.Material.dEdx(beta,projectedData(9,:),projectedData(8,:));
            % positions of samples on z-theta*R graph
            z = projectedData(3,:);
            theta = sign(projectedData(2,:)).*acos(projectedData(1,:)./sqrt(projectedData(1,:).^2+projectedData(2,:).^2));
            y = theta*Consts.TubeRadius;
            % store zysamples on the GPU to speed it up for larger samples
            zySamples = [z;y]';


            yGrid = equiprobableGrid(y,options.MaxThetaPoints);

            if options.echo
                % begin performance timer
                fprintf("Performing kernel energy density estimation\n");
                tic
            end
            % perform 3-mean clustering on log vals of z axis
            logZ = log10(abs(z));
            % TODO make k and optional parameter
            k = options.kmeans;
            if k==1
                clusterIdx = ones(1,length(z));
            else
                clusterIdx = kmeans(logZ',k)';
            end

            cellSamples = cell(1,k);
            cellWeights = cell(1,k);
            cellGrids = cell(1,k);
            maximumIdxInCluster = zeros(1,k);
            maximumValInCluster = zeros(1,k);
            maximumCombinedEnergies = zeros(1,k);
            for j = 1:k
                thisZGrid = equiprobableGrid(z(clusterIdx==j),options.MaxZPoints);
                cellSamples(j) = {zySamples(clusterIdx==j,:)};
                cellWeights(j) = {dEdx(clusterIdx==j)};
                [Zj,Yj] = meshgrid(thisZGrid,yGrid);
                cellGrids(j) = {[Zj(:) Yj(:)]};
            end
            % perform kde on each cluster:
            parfor i = 1:k
                % find total energy of the cluster
                clusterdEdx = cellWeights{i};
                clusterEnergy = sum(clusterdEdx);
                % make weights from relative energy to cluster total
                clusterWeights = clusterdEdx./clusterEnergy';
                % apply weighted kde to global grid

                clusterPDF = clusterEnergy.*ksdensity(cellSamples{i},cellGrids{i},"weights",clusterWeights);
                % get cluster's contribution to total energy on grid
                [maxVal,I] = max(clusterPDF);
                maximumValInCluster(i) = maxVal;
                maximumIdxInCluster(i) = I;
            end
            if options.echo
                tElapsed = toc;
                seconds = rem(tElapsed,60);
                minutes = (tElapsed-seconds)./60;
                fprintf("Energy deposition calculated in %d min %.5f\n",minutes,seconds)
            end

            % evaluate and sum the kdes at each maximum
            for i = 1:k
                currTotal = 0;
                currGrid = cellGrids{i};
                maximumLocation = currGrid(maximumIdxInCluster(i),:);
                for j = 1:k
                    if i==j
                        currTotal = currTotal + maximumValInCluster(i);
                    else
                        currTotal = currTotal + sum(cellWeights{j}).*ksdensity(cellSamples{j},maximumLocation,"Weights",cellWeights{j}./sum(cellWeights{j}));
                    end
                end
                maximumCombinedEnergies(i) = currTotal;
            end
            E = max(maximumCombinedEnergies)*1e9*Consts.e;

        end

        function [primaryRange,secondaryRange] = getRange(obj,material,options)
            arguments(Input)
                obj(1,1) MCResult
                material(1,1) Material = "W"
                options.include(1,1) string {mustBeMember(options.include,["primaries","secondaries","all"])} = "all"
                options.method(1,1) string {mustBeMember(options.method,["True","Perpendicular"])} = "Perpendicular"
            end
            [primaryRange,secondaryRange] = obj.getCSDARange(material,include=options.include);
            % mean angle follows rayleigh distribution with theta_0

            
            
            Z = material.Z;
            A = material.A;
            wt = material.wt;
            X0 = material.X0;
            rho = material.rho;
            primaryBeta = sqrt(1-1./(obj.projectedPrimaries(4,:)./obj.projectedPrimaries(9,:)).^2);
            primaryMomentum = obj.projectedPrimaries(4,:).*primaryBeta;
            primaryCharge = obj.projectedPrimaries(8,:);
            parfor i = 1:length(primaryRange)
                primaryTheta = sqrt(pi/2) * MCS.getTheta0(primaryMomentum(i),primaryBeta(i),primaryCharge(i),Z,A,X0,rho,wt,primaryRange(i));
                primaryRange(i) = 2*primaryRange(i)./(1+1/cos(primaryTheta));
            end

            secondaryBeta = sqrt(1-1./(obj.projectedSecondaries(4,:)./obj.projectedSecondaries(9,:)).^2);
            secondaryMomentum = obj.projectedSecondaries(4,:).*secondaryBeta;
            secondaryCharge = obj.projectedSecondaries(8,:);
            parfor i = 1:length(secondaryRange)
                secondaryTheta = sym(sqrt(pi/2) * MCS.getTheta0(secondaryMomentum(i),secondaryBeta(i),secondaryCharge(i),Z,A,X0,rho,wt,secondaryRange(i)));
                secondaryRange(i) = double(2*sym(secondaryRange(i))./(1+1/cos(secondaryTheta)));
            end

            if options.method=="Perpendicular"
                if options.include == "primaries"||options.include=="all"
                    primaryDirections = sym(obj.projectedPrimaries(5:7,:));
                    primarySines = double(sin(acos(primaryDirections(3,:)./sqrt(primaryDirections(3,:).^2+primaryDirections(2,:).^2+primaryDirections(1,:).^2))));
                    primaryRange = primaryRange.*abs(primarySines);
                end
                if options.include=="secondaries"||options.include=="all"
                    secondaryDirections = sym(obj.projectedSecondaries(5:7,:));
                    secondarySines = double(sin(acos(secondaryDirections(3,:)./sqrt(secondaryDirections(3,:).^2+secondaryDirections(2,:).^2+secondaryDirections(1,:).^2))));
                    secondaryRange = secondaryRange.*secondarySines;
                end
            end
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