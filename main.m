%% Define geometry
nSamples = 1e5;
sigma = 220e-6;
thickness = 1e-2;
mu = 0.002;

haloLimit = 4.8*sigma;
offset = haloLimit+mu;
runtimeData = Consts.data;
runtimeData.default;

% SAMPLES
% Samples are stored as column vectors of length 8
% entries 1-3 store 3d position
% entries 4-7 store 3d four momentum  in natural units GeV/c
% entry 8 stores particle charge
% entry 9 stores rest mass in natural units GeV/c^2 for numerical accuracy
% geometry = SymRectangles(0.01,0.01,2*offset,thickness);
geometry = OffsetRectangle(0.01,0.01,offset,thickness);
geometry.xRotation = 0;
[pdf,Ihalo,expectedProtons] = constantHalo(0,sigma,1e-5,geometry.xLimits,geometry.yLimits);
% f = normalDistribution(0,220);
% expectedProtons = 1/(sigma^2*2*pi)*exp(-0.5*offset^2./sigma^2)*Consts.data.Nb/Consts.data.bunchSpacing*30e-6*0.01;
% geometry.pdf = pdf;
% scale factor = expected incident protons per second/simulated protons
scaleFactor = expectedProtons/nSamples;

material = Material("Cu");
% make some samples (pdf is automatically set to uniform)
samples = geometry.generate(nSamples);

%%% Prepare simulation
% make the input object
mcinput = MCInput(geometry,material,samples=samples);
% make simulation options object
%%
simopts = SimulationOptions(exclude="all",threads=8,savesPositions=false,recursionLimit=2);
simopts.include("Ionisation","Nuclear")
% mcinput.view
% hold on
% plot(0,0,"k.")
% plot(mu,0,"m*")
%% Run simulation
mcout = mcsimulate(mcinput,simopts);
%% Inspect results
mcout.viewIntersections("scalefactor",scaleFactor,"MaxThetaPoints",30);