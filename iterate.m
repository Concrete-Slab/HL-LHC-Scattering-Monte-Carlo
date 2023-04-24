%% SETUP
sigma = 7.1e-6;
mu = 300e-6;

haloLimit = 4.8*sigma;
offset = haloLimit+mu;
screenSideLength = 0.01;
xLimits = [offset;offset + screenSideLength];
yLimits = [-screenSideLength/2;screenSideLength/2];
runtimeData = Consts.data;
runtimeData.default;
M = Material("Chromox");
[pdf,Ihalo,expectedProtons] = constantHalo(mu,sigma,1e-5,xLimits,yLimits);
simopts = SimulationOptions(exclude="all",threads=8,savesPositions=false,recursionLimit=2);
simopts.include("MCS","Bethe","Ionisation","Rutherford")
%% RUN SOLVER
f = @(t) runSim(M,simopts,t,screenSideLength,offset,pdf,expectedProtons);


%[tGuess,Eguess] = solveBisection(f,5e-5,1e-2,1e-6,20);
str = ["YAGCe"];%,"NaITl","BaF2","CeBr3","Al","Si","Graphite","Cu","Ti"];
mats = repmat(Material,1,length(str));
for i = 1:length(mats)
    mats(i) = Material(str(i));
end

maxThickness = zeros(1,length(mats));
deposition = zeros(1,length(mats));
for i = 1:length(mats)
    fprintf("\nCurrent material: %s (%d/%d)\n",str(i),i,length(mats))
    tMaterial = tic;
    f = @(t) runSim(mats(i),simopts,t,screenSideLength,offset,pdf,expectedProtons);
    [Tguess,Eguess] = solveBisection(f,1e-4,2e-4,1e-6,20);
    maxThickness(i) = Tguess;
    deposition(i) = Eguess;
    elapsedTime = toc(tMaterial);
    seconds = rem(elapsedTime,60);
    minutes = (elapsedTime-seconds)./60;
    fprintf("Finished calculation for %s, took %d mins %.5f seconds\n",str(i),minutes,seconds)
end


%%
load chromomxThickness.mat
mats2 = [mats Material("Chromox")];
str2 = [str "Chromox"];
deposition2 = [deposition Eguess];
maxThickness2 = [maxThickness tGuess];
[maxThicknessSorted, I] = sort(maxThickness2);
depositionSorted = deposition2(I);
matsSorted = mats2(I);
strSorted = str2(I);
OTR = ismember(strSorted,["Al","Ti","Si","Graphite","Cu"]);

hold on
yyaxis left
b = bar(1:length(strSorted),maxThicknessSorted*1e6);
ax = gca;
b.FaceColor='flat';
b.CData(OTR,:) = repmat([0,0,1],length(find(OTR)),1);
b.CData(~OTR,:) = repmat([0,0.4,0.6],length(find(~OTR)),1);
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



function val = runSim(material,simopts,thickness,sideLength,offset,pdf,expectedProtons)
    log10nMax = 7;
    log10nMin = 4;
    log10tMax = 6;
    log10tMin = 2;
    nSamples = ceil(10^((log10nMax-log10nMin)./(log10tMax-log10tMin) * (-log10(thickness)-log10tMin)+log10nMin));
%     nSamples = ceil(10^(-2/3*log10(thickness)+8/3));
    scaleFactor = expectedProtons/nSamples;
    geometry = OffsetRectangle(sideLength,sideLength,offset,thickness,0);
    geometry.pdf = pdf;
    samples = geometry.generate(nSamples,"echo","off");
    in = MCInput(geometry,material,samples=samples);
    out = mcsimulate(in,simopts,"echo","off");
    [~,~,E] = out.energyDistribution("MaxThetaPoints",20,"echo","off");
    E = E.*scaleFactor.*1e-3;
    val = max(max(E))-4.5;
end