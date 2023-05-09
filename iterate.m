%% SETUP
sigma = 80e-6;
mu = 0.002;
nSamples = 1e6;
t10 = 1e-2;
t00 = 1e-5;
maxIterations = 25;
abstol = 2.5e-6;
screenSideLength = 0.01;
simopts = SimulationOptions(exclude="all",threads=8,savesPositions=false,recursionLimit=0);
simopts.include("Nuclear","MCS","Ionisation","Bethe","Rutherford")
kmeans = 2;
scaleFactor = Consts.data.Nb./Consts.data.bunchSpacing * 10^-5/(2*pi*sigma^2) * screenSideLength^2 * 1e-3/nSamples;
%% RUN SOLVER
str = ["YAPCe","LSOCe"];
mats = repmat(Material,1,length(str));
for i = 1:length(mats)
    mats(i) = Material(str(i));
end
%%
maxThickness = zeros(1,length(mats));
deposition = zeros(1,length(mats));
for i = 1:length(mats)

    fprintf("\nCurrent material: %s (%d/%d)\n",str(i),i,length(mats))
    tMaterial = tic;
    t1 = t10;
    ft1 = 1;
    t0 = t00;
    ft0 = -1;
    if ft0*ft1>0
        error("Both guesses have same sign")
    end
    j=0;
    
    while j<=maxIterations
        j = j+1;
        t2 = 10^((log10(t1)+log10(t0))/2);
        if abs(t0-t2)<abstol
            break
        end
        tStart = tic;
            g = OffsetRectangle(screenSideLength,screenSideLength,4.8*sigma+mu,t2);
            out = mcsimulate(MCInput(g,mats(i),"n",nSamples,"echo","off"),simopts,echo="off");
            ft2 = scaleFactor*out.maxEnergy("kmeans",kmeans,"MaxThetaPoints",20,"echo","off") - 4;
        tStop = toc(tStart);
        seconds = rem(tStop,60);
        minutes = (tStop-seconds)./60;
        fprintf(" - Guess %d is %.1f microns. Sim time was %d min, %.5fs\n",j,t2*1e6,minutes,seconds)
        if ft2*ft0>0
            % new value lies between t2 and t1
            t0 = t2;
            ft0 = ft2;
        else
            t1 = t2;
        end
    end
    



    elapsedTime = toc(tMaterial);
    seconds = rem(elapsedTime,60);
    minutes = (elapsedTime-seconds)./60;
    fprintf("Finished calculation for %s, took %d mins %.5f seconds. Thickness is %.1f microns\n",str(i),minutes,seconds,t2*1e6)
    maxThickness(i) = t2;
    deposition(i) = ft2 + 4;
end
% save MaxThicknesses mats str maxThickness deposition
save ExtraThicknesses mats str maxThickness deposition
%%
[maxThicknessSorted, I] = sort(maxThickness);
depositionSorted = deposition(I);
matsSorted = mats(I);
strSorted = str(I);
OTR = ismember(strSorted,["Al","Ti","Si","Graphite","Cu"]);

hold on
yyaxis left
b = bar(1:length(strSorted),maxThicknessSorted*1e6);
ax = gca;
b.FaceColor='flat';
b.CData(OTR,:) = repmat([0,0,1],length(find(OTR)),1);
b.CData(~OTR,:) = repmat([0,1,0],length(find(~OTR)),1);
ylabel("Maximum thickness, \mu{m}")
yyaxis right
merit = zeros(1,length(matsSorted));
merit2 = merit;
for j = 1:length(merit)
    m = matsSorted(j);
    merit(j) = 1./(sum(m.wt.*m.Z.*(m.Z+1)./m.A).*m.rho);
    merit2(j) = merit(j)./m.braggI;
end
merit2 = merit2.*(merit(end)/merit2(end));
plot(1:length(strSorted),merit)
%plot(1:length(strSorted),merit2)
ylabel("Scattering merit M")
set(ax,"XTickLabel",strSorted)
set(ax,"Xtick",1:length(strSorted))
xlabel("Material")



% function val = runSim(material,simopts,thickness,sideLength,offset,pdf,expectedProtons)
%     log10nMax = 7;
%     log10nMin = 4;
%     log10tMax = 6;
%     log10tMin = 2;
%     nSamples = ceil(10^((log10nMax-log10nMin)./(log10tMax-log10tMin) * (-log10(thickness)-log10tMin)+log10nMin));
% %     nSamples = ceil(10^(-2/3*log10(thickness)+8/3));
%     scaleFactor = expectedProtons/nSamples;
%     geometry = OffsetRectangle(sideLength,sideLength,offset,thickness,0);
%     geometry.pdf = pdf;
%     samples = geometry.generate(nSamples,"echo","off");
%     in = MCInput(geometry,material,samples=samples);
%     out = mcsimulate(in,simopts,"echo","off");
%     [~,~,E] = out.energyDistribution("MaxThetaPoints",20,"echo","off");
%     E = E.*scaleFactor.*1e-3;
%     val = max(max(E))-4.5;
% end