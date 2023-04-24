function out = mcsimulate(mcinput,simopt,options)
%SIM Summary of this function goes here
%   Detailed explanation goes here
arguments(Input)
    mcinput(1,1) MCInput
    simopt(1,1) SimulationOptions = SimulationOptions
    options.echo(1,1) string {mustBeMember(options.echo,["on","off"])} = "on";
end
arguments(Output)
    out(1,1) MCResult
end
options.echo = ismember(options.echo,"on");
% read input
samplesIn = mcinput.samples;
nSamples = mcinput.nSamples;
geom = mcinput.geometry;

% configure options
nThreads = simopt.threads;
[hard,soft] = simopt.configure(mcinput);
savesPositions = simopt.savesPositions;
maxRecursion = simopt.recursionLimit;
maxSecondariesPerRecursion = simopt.estimatedSecondaries(mcinput);


% preallocate outputs
primariesOut = zeros(9,nSamples);
secondariesOut = repmat(SecondarySamples,1,nSamples);
if savesPositions
    primaryHistories = repmat(ParticleHistory,1,nSamples);
else
    primaryHistories = ParticleHistory.empty;
end

if nThreads>1
    % create a parallel pool if not already done, with the desired thread num
    if isempty(gcp('nocreate'))
        pool = parpool(nThreads);
    else
        pool = gcp('nocreate');
    end
    threadsUsed = pool.NumWorkers;
else
    threadsUsed = 1;
end

if options.echo
    % prepare output material
    fprintf("Beginning simulation on %d threads, with the following parameters:\n",threadsUsed)
    if simopt.savesPositions
        fprintf("Particle intermediate position saving - True\n")
    else
        fprintf("Particle intermediate position saving - False\n")
    end
    if isempty(simopt.included)
        strout = "None";
    else
        strout = simopt.included;
    end
    fprintf("Models included - %s\n",strout.join(", "))
    fprintf("Maximum secondary simulation recursion depth - %d\n",maxRecursion)
    % begin performance timer
    tic
end


if threadsUsed>1
    % begin main loop in parallel
    parfor i=1:nSamples
        [primary,secondaries,particleHistory,secondaryHistory] = simulateOne(samplesIn(:,i),hard,soft,geom,savesPositions,0,maxRecursion,maxSecondariesPerRecursion);
        primariesOut(:,i) = primary;
        if savesPositions
            primaryHistories(i) = particleHistory;
        end
        secondariesOut(i) = SecondarySamples(secondaries,secondaryHistory);
    end
else
    % main loop is regular for loop if only one thread
    for i=1:nSamples
        [primary,secondaries,particleHistory,secondaryHistory] = simulateOne(samplesIn(:,i),hard,soft,geom,savesPositions,0,maxRecursion,maxSecondariesPerRecursion);
        primariesOut(:,i) = primary;
        if savesPositions
            primaryHistories(i) = particleHistory;
        end
        secondariesOut(i) = SecondarySamples(secondaries,secondaryHistory);
    end
end
if options.echo
    % save elapsed simulation time
    elapsedSimulation = toc;
    secondsSim = rem(elapsedSimulation,60);
    minutesSim = (elapsedSimulation-secondsSim)./60;
    % display simulation time and info
    fprintf("Simulation completed in %d min %.5f s\n",minutesSim,secondsSim)
    tic
end
% configure output object and data
out = MCResult(primariesOut,secondariesOut,mcinput,primaryHistories);
% display time taken to configure output
if options.echo
    elapsedOutput = toc;
    secondsOut = rem(elapsedOutput,60);
    minutesOut = (elapsedOutput-secondsOut)./60;
    fprintf("Output object configured in %d min %.5f s\n",minutesOut,secondsOut);
end
end


% function [primary,secondaries,primaryHistory,secondaryHistory] = simulateOne(initialSample,hard,soft,geometry,savesPositions,recursionDepth,maxRecursion,maxSecondariesPerRecursion)
% arguments(Input)
%     initialSample(9,1) double
%     hard(1,:) HardProcess
%     soft(1,:) SoftProcess
%     geometry(1,1) Geometry
%     savesPositions(1,1) logical
%     recursionDepth(1,1) double {mustBeInteger,mustBeNonnegative}
%     maxRecursion(1,1) double {mustBeInteger,mustBeNonnegative}
%     maxSecondariesPerRecursion(1,1) double {mustBeInteger,mustBeNonnegative}
% end
% arguments(Output)
%     primary(9,1) double
%     secondaries(9,:) double
%     primaryHistory(1,:) ParticleHistory {mustBeScalarOrEmpty}
%     secondaryHistory(1,:) ParticleHistory
% end
% % initialise primary particle's properties
% pfm = initialSample(4:7);
% penergy = pfm(1);
% pmomentum = pfm(2:4);
% ppos = initialSample(1:3);
% pcharge = initialSample(8);
% pmass = initialSample(9);
% 
% % set up a particle handle so that models can all listen to updates
% particle = ParticleHandle(penergy,pcharge,pmass);
% 
% if savesPositions
%     % set up a particle history object to track the primary particle
%     % estimate a max number of interactions to be 100
%     primaryHistory = ParticleHistory(100,ppos);
% else
%     primaryHistory = ParticleHistory.empty;
% end
% nhard = length(hard);
% nsoft = length(soft);
% for i = 1:nhard
%     hard(i).particle = particle;
% end
% for j = 1:nsoft
%     soft(j).particle = particle;
% end
% 
% % preallocate the arrays to store secondaries
% if maxSecondariesPerRecursion>0
%     sizeAtCurrentDepth = maxSecondariesPerRecursion^(maxRecursion-recursionDepth);
%     secondaries = nan(9,sizeAtCurrentDepth);
%     if savesPositions
%         secondaryHistory = repmat(ParticleHistory,1,sizeAtCurrentDepth);
%     else
%         secondaryHistory = ParticleHistory.empty(1,0);
%     end
% else
%     secondaries = double.empty(9,0);
%     secondaryHistory = ParticleHistory.empty(1,0);
% end
% secondaryIndex = 1;
% 
% % begin loop
% propagate = true;
% while propagate
% 
%     withinBounds = false;
%     direction = pmomentum./norm(pmomentum);
%     if nhard>0
%         % update delta x, then calculate new angles etc
% 
%         % First, recalculate the total inverse mean free path
%         M0 = zeros(1,nhard);    % individual inverse mean paths
%         cmf = zeros(1,nhard);   % cumulative mass function of M0(i)
%         for i = 1:nhard
%             M0(i) = hard(i).M0; % calculate M0 from hard process
%             if i>1              % update cmf if there is more than one hard process
%                 cmf(i) = cmf(i-1) + M0(i);
%             else
%                 cmf(i) = M0(i);
%             end
%         end
%         invlambda = sum(M0);
%         if invlambda ~=0
%             cmf = cmf./invlambda;
%         else
%             % all processes have M0 of zero! (dummy list)
%             % equal chance to choose each ( although 
%             cmf = (1:nhard)./nhard;
%         end
% 
%         % Draw a random deltaX from exponential distribution
%         zDx = rand;
%         deltaX = -log(zDx)./invlambda;
% 
%         % check if moving by deltaX puts particle outside the geometry
%         projectedPosition = ppos+deltaX*direction;
%         withinBounds = geometry.isInside(projectedPosition);
% 
%         if withinBounds % particle will undergo a hard interaction
% 
%             % move particle into position
%             ppos = projectedPosition;
% 
%             % randomly select process from cmf
%             zP = rand;
%             possibleProcesses = hard(cmf>=zP);
%             hProcess = possibleProcesses(end);
% 
%             % perform the interaction
%             % the interact command should automatically update the particle
%             % handle, so every process should register the update
%             [dE,dAnglesHard,newSecondaryInfo] = hProcess.interact();
%             % rotation matrix to convert between coordinate rotation of
%             % interaction (particle initially along z) and global frame
%             R = [[1;0;0] [0;1;0] direction];
%             newSecondaryInfo(2:4) = R*newSecondaryInfo(2:4);
%             % the hard interaction may create a secondary particle!
%             % this must be simulated through recursion of this function
%             % the conditions for further simulation of the secondary are:
%             %   the hard process must have created a secondary
%             %   AND simulating will not increase recursion past the maximum
% 
% 
%             if recursionDepth<maxRecursion && size(newSecondaryInfo,2)==1 && secondaryIndex<sizeAtCurrentDepth
%                 % recurse with the same parameters, except:
%                 %   primary particle is now the secondary from the process
%                 %   recursion depth is increased by 1
%                 newSecondaryParticle = [ppos;newSecondaryInfo];
%                 [initSecondary,additionalSecondaries,initHistory,additionalHistories] = simulateOne(newSecondaryParticle,hard,soft,geometry,savesPositions,recursionDepth+1,maxRecursion,maxSecondariesPerRecursion);
%                 % this line is not efficient - maybe estimate sizes for a
%                 % given maximum depth?
%                 nNewSecondaries = length(additionalHistories) + 1;
%                 nextSecondaryIndex = secondaryIndex + nNewSecondaries - 1;
%                 if nextSecondaryIndex<sizeAtCurrentDepth
%                     secondaries(:,secondaryIndex:nextSecondaryIndex) = [initSecondary,additionalSecondaries];
%                     if ~savesPositions
%                         % secondaries need at least 2 history points so
%                         % they can be followed. initHistory will be empty
%                         % because it is treated as the primary particle
%                         initHistory = ParticleHistory(2);
%                         initHistory = initHistory.write([ppos,initSecondary(1:3)]);
%                     end
%                     secondaryHistory(:,secondaryIndex:nextSecondaryIndex) = [initHistory,additionalHistories];
%                 end
%                 secondaryIndex = nextSecondaryIndex + 1;
% 
%             end
% 
%             % save the new momentum/energy
%             penergy = particle.energy;
% 
%             % assume radial symmetry in scattering
%             % the inclination angle theta is determined by interact()
%             % the polar angle phi is randomly sampled between 0 and 2*pi
% %             dPhi = 2*pi*rand;
% 
% 
% 
%             % transform new direction in particle frame to global frame
%             direction = R*[dAnglesHard(2)*dAnglesHard(3) ; dAnglesHard(1)*dAnglesHard(3) ; dAnglesHard(4)];
%             direction = direction./norm(direction);
%             pmomentum = particle.momentum*direction;
% 
% %             if direction(2)==0&&direction(1)==0
% %                 initialPhi=0;
% %             else
% %                 initialPhi = sign(pmomentum(2))*acos(pmomentum(1)./sqrt(pmomentum(1)^2+pmomentum(2)^2));
% %             end
% %             phi = initialPhi+dPhi;
% % 
% %             % use polar relations to switch back to cartesian momentum
% %             pmomentum = momentumMagnitude.*[cos(phi)*sin(theta) ; sin(phi)*sin(theta) ; cos(theta)];
%             pfm = [penergy;pmomentum];
%         end
%     end
% 
%     if (~withinBounds||nhard==0)
%         % particle traverses remaining distance unaffected by hard process
%         propagate = false;
% %         if ~geometry.isInside(ppos)
% %             ppos
% %         end
%         projectedPosition = geometry.edgeIntersect(ppos,direction);
%         deltaX = norm(projectedPosition-ppos);
%         % update particle position only (energy,momentum are same)
%         ppos = projectedPosition;
%     end
% 
% 
%     % apply soft processes after hard
%     if nsoft>0
%         % apply soft processes using delta x
%         dE = 0;
%         dAnglesSoft = 0;
%         dPos = 0;
%         for i = 1:nsoft
%             [dEi,dAnglesi,dPosi] = soft(i).update(deltaX);
%             % assume all soft processes are independent, so add
%             % contributions
%             dE = dE+dEi;
%             dAnglesSoft = dAnglesSoft + dAnglesi;
%             dPos = dPos + dPosi;
%         end
% 
%         % make this more efficient?
%         
%         % save the new position/momentum/energy
% 
%         % dPos is only in transverse direction, hence select 1:2
%         ppos(1:2) = ppos(1:2) + dPos;
%         penergy = particle.energy;
% 
%         % rotation matrix method
%         R = [[1;0;0] [0;1;0] direction];
%         direction = R*[cos(dAnglesSoft(1))*sin(dAnglesSoft(2)) ; sin(dAnglesSoft(1))*sin(dAnglesSoft(2)) ; cos(dAnglesSoft(2))];
%         direction = direction./norm(direction);
%         pmomentum = particle.momentum*direction;
%         pfm = [penergy;pmomentum];
%     else
%         % do nothing
%     end
%     % save the position of the particle
%     if savesPositions
%         primaryHistory = primaryHistory.write(ppos);
%     end
% end
% 
% % clean up memory from the particle handle
% particle.delete();
% 
% if savesPositions
%     % remove zeros in particle history to reduce memory transfer size
%     primaryHistory = primaryHistory.memsave();
% end
% primary = [ppos;pfm;pcharge;pmass];
% if any(~isfinite(primary))
%     primary
% end
% secondaries = secondaries(:,all(~isnan(secondaries)));
% secondaryHistory = secondaryHistory(all(~isnan(secondaries)));
% 
% end
