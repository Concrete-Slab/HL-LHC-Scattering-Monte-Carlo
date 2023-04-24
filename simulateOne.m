function [primary,secondaries,primaryHistory,secondaryHistory] = simulateOne(initialSample,hard,soft,geometry,savesPositions,recursionDepth,maxRecursion,maxSecondariesPerRecursion)
arguments(Input)
    initialSample(9,1) double
    hard(1,:) HardProcess
    soft(1,:) SoftProcess
    geometry(1,1) Geometry
    savesPositions(1,1) logical
    recursionDepth(1,1) double {mustBeInteger,mustBeNonnegative}
    maxRecursion(1,1) double {mustBeInteger,mustBeNonnegative}
    maxSecondariesPerRecursion(1,1) double {mustBeInteger,mustBeNonnegative}
end
arguments(Output)
    primary(9,1) double
    secondaries(9,:) double
    primaryHistory(1,:) ParticleHistory {mustBeScalarOrEmpty}
    secondaryHistory(1,:) ParticleHistory
end
% initialise primary particle's properties
mke = Consts.MinimumKineticEnergy;
pfm = initialSample(4:7);
penergy = pfm(1);
pmomentum = pfm(2:4);
ppos = initialSample(1:3);
pcharge = initialSample(8);
pmass = initialSample(9);
minimumMaximumStepLength = geometry.thickness + eps;
% set up a particle handle so that models can all listen to updates
particle = ParticleHandle(penergy,pcharge,pmass);

if savesPositions
    % set up a particle history object to track the primary particle
    % estimate a max number of interactions to be 100
    primaryHistory = ParticleHistory(100,ppos);
else
    primaryHistory = ParticleHistory.empty;
end
nhard = length(hard);
nsoft = length(soft);
for i = 1:nhard
    hard(i).particle = particle;
end
for i = 1:nsoft
    soft(i).particle = particle;
end

% preallocate the arrays to store secondaries
if maxSecondariesPerRecursion>0
    sizeAtCurrentDepth = maxSecondariesPerRecursion^(maxRecursion-recursionDepth);
    secondaries = nan(9,sizeAtCurrentDepth);
    if savesPositions
        secondaryHistory = repmat(ParticleHistory,1,sizeAtCurrentDepth);
    else
        secondaryHistory = repmat(ParticleHistory(2),1,sizeAtCurrentDepth);
    end
else
    secondaries = double.empty(9,0);
    secondaryHistory = ParticleHistory.empty(1,0);
end
secondaryIndex = 1;



% begin loop
propagate = true;
direction = pmomentum./norm(pmomentum);
% GENERATE TRANSFORM TO LOCAL COORDINATES
T = localTransform(direction);
hardIndices = 0:nhard;
while propagate
    dEdx = 0;
    for i = 1:nsoft
        dEdx = dEdx + soft(i).dEdx;
    end
    % dont allow particle to lose more than 2% of energy.
    % if this step size is too small, then use 1e-7
    maxStep = min(minimumMaximumStepLength,0.02*particle.energy./dEdx);
    hardInteraction = false;
    
    
    if nhard>0
        % initially sample a delX
        M0 = zeros(1,nhard);
        cmf = zeros(1,nhard+1);
        % individual inverse mean free paths
        M0(1) = hard(1).M0;
        cmf(2) = M0(1);
        for i = 2:nhard
            M0(i) = hard(i).M0;
            cmf(i+1) = cmf(i)+M0(i);
        end
        inverselambda = sum(M0);
        % cumulative mass function for sampling from
        if inverselambda~=0
            cmf = cmf./inverselambda;
        else
            cmf = (1:nhard)./nhard;
        end

        % sample a random delX from the exponential distribution
        % Draw a random deltaX from exponential distribution
        zDx = rand;
        deltaX = -log(zDx)/inverselambda;
        if deltaX<maxStep
            hardInteraction = true;
        else
            deltaX=maxStep;
        end
    else
        deltaX = maxStep;
    end
    
    % prepare to move by delX
    projectedPosition = ppos + deltaX*direction;
    % if delX is too large, then bring to edge, cease propagation and do
    % not apply a hard interaction, as the particle never reaches
    % interaction point
    if ~geometry.isInside(projectedPosition)
        projectedPosition = geometry.edgeIntersect(ppos,direction);
        
        deltaX = norm(projectedPosition-ppos);
        propagate = false;
        hardInteraction = false;
    end

    % first apply soft processes
    dAngleSoft = [0;0];
    dPos = [0;0];
    dE = 0;
    for i = 1:nsoft
        [dEi,dAngleSofti,dPosi] = soft(i).update(deltaX);
        dE = dE+dEi;
        dAngleSoft = dAngleSoft + dAngleSofti;
        dPos = dPos + dPosi;
    end

    if particle.kineticEnergy+dE<mke
        % hard process now cannot happen because particle never reaches
        % interaction point
        hardInteraction=false;
        propagate = false;
        particle.momentum = 0;
    else
        particle.energy = particle.energy-dE;
    end

    % move the particle - intermediate diagram
    if ~all(dPos==0) % only perform below code if position changes softly
        proposedPosition = projectedPosition + T*[dPos;0];
        if ~geometry.isInside(proposedPosition)
            % proposed new position is outside of geometry
            % we need to interpolate, position, stop the propagation and
            % prevent a hard interaction from taking place
            propagate = false;
            hardInteraction = false;
            % position + x/2
            intermediatePosition = ppos + 0.5*(projectedPosition-ppos);
            propagationDirection = proposedPosition-intermediatePosition;
            finalPosition = geometry.edgeIntersect(intermediatePosition,propagationDirection./norm(propagationDirection));
            ppos = finalPosition;
        else
            ppos = proposedPosition;
        end
    else % otherwise just move to the projected position
        ppos = projectedPosition;
    end
%     ppos = projectedPosition + T*[dPos;0];
%     if any(~isreal(ppos))
%         ppos
%     end
%     if ~geometry.isInside(ppos)
%         propagate = false;
%         hardInteraction=false;
%     end
    % change its direction
    newDirection = T*[cos(dAngleSoft(1))*sin(dAngleSoft(2));sin(dAngleSoft(1))*sin(dAngleSoft(2));cos(dAngleSoft(2))];
    direction = newDirection./norm(newDirection);
    T = localTransform(direction);
    
    if hardInteraction
        
        if nhard>1
            % randomly select hard process by probability mass
            zP = rand;
            [~,Ia,~] = unique(cmf,"first");
            hi = hardIndices(Ia);
            iSampled = interp1(cmf(Ia),hi,zP);
            iChosen = min(hi(hi>iSampled));
            hp = hard(iChosen);
        else
            hp = hard(1);
        end
        
        % undergo hard interaction
        [dE,dAngleHard,newSecondaryInfo] = hp.interact;
        % check if the particle is still moving enough
        if particle.kineticEnergy+dE<mke
            particle.momentum = 0;
            propagate = false;
        else
            particle.energy = particle.energy + dE;
        end
        newDirection = T*[dAngleHard(2)*dAngleHard(3);dAngleHard(1)*dAngleHard(3);dAngleHard(4)];
        direction = newDirection./norm(newDirection);
        T = localTransform(direction);
        if ~isempty(newSecondaryInfo) && recursionDepth<maxRecursion && secondaryIndex<sizeAtCurrentDepth
            % follow the secondary with recursion
            newSecondaryInfo(2:4) = T*newSecondaryInfo(2:4);
            [initSecondary,additionalSecondaries,initHistory,additionalHistories] = simulateOne([ppos;newSecondaryInfo],hard,soft,geometry,savesPositions,recursionDepth+1,maxRecursion,maxSecondariesPerRecursion);
            
            % do we have enough room to store all the new secondaries?
            nSecondaries = length(additionalHistories)+1;
            nextSecondaryIndex = secondaryIndex+nSecondaries-1;
            if nextSecondaryIndex<=sizeAtCurrentDepth
                secondaries(:,secondaryIndex:nextSecondaryIndex) = [initSecondary,additionalSecondaries];
                if ~savesPositions
                    % secondaries need at least 2 history points so
                    % they can be followed. initHistory will be empty
                    % because it is treated as the primary particle
                    initHistory = ParticleHistory(2);
                    initHistory = initHistory.write([ppos,initSecondary(1:3)]);
                end
                secondaryHistory(:,secondaryIndex:nextSecondaryIndex) = [initHistory,additionalHistories];
            end
            secondaryIndex = nextSecondaryIndex+1;
        end
    end

    % save the position of the particle
    if savesPositions
        primaryHistory = primaryHistory.write(ppos);
    end
end
pfm = [particle.energy;particle.momentum.*direction];
% clean up memory from the particle handle
particle.delete();

if savesPositions
    % remove zeros in particle history to reduce memory transfer size
    primaryHistory = primaryHistory.memsave();
end

primary = [ppos;pfm;pcharge;pmass];
% if ~geometry.isInside(ppos)
%     primary
% end
actualSecondaries = all(~isnan(secondaries));
secondaries = secondaries(:,actualSecondaries);
secondaryHistory = secondaryHistory(actualSecondaries);

end