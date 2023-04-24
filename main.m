%% Define geometry
nSamples = 1e1;
sigma = 7.1e-6;
thickness = 1e-2;
mu = 300e-6;

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

geometry.pdf = pdf;
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
simopts = SimulationOptions(exclude="all",threads=8,savesPositions=true,recursionLimit=2);
simopts.include("Bethe","Ionisation","Rutherford")
% mcinput.view
% hold on
% plot(0,0,"k.")
% plot(mu,0,"m*")
%% Run simulation
mcout = mcsimulate(mcinput,simopts);
%% Inspect results
mcout.viewIntersections("scalefactor",scaleFactor,"MaxThetaPoints",30);