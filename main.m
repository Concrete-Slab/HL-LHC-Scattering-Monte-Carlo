Consts.data.default
% IMPORTANT - add subdirectories to MATLAB search path with below command
addpath(genpath(pwd))
%% Define problem geometry
nSamples = 1e2; % number of protons that will be simulated

% HL-LHC beam properties
sigma = 220e-6; % beam rms width, metres
mu = 0.002; % distance between beam centre and tube centre, metres


% screen geometry properties
thickness = 1e-3; % thickness of screen, metres
sideLength = 1e-2; % side lengths of square screen, metres


% make the geometry object
offset = 4.8*sigma+mu; % max. metres from tube centre where halo begins
geometry = OffsetRectangle(sideLength,sideLength,offset,thickness);


% assuming constant halo and square screen:
% - find the number of protons passing through the screen in the HL-LHC
% - find the factor by which this is greater than the number of samples
Ihalo = getHaloIntensity(sigma);
expectedProtons = Ihalo*geometry.area;
scaleFactor = expectedProtons/nSamples;


% choose a material for the screen
% supported string names for materials are given by Materials.allNames
material = Material("LSOCe");

%% Prepare simulation objects
% create an input object, generate samples
mcinput = MCInput(geometry,material,n=nSamples);

% configure options for the simulation
% - set the # of threads: for performance use up to 8, for debugging use 1
% - choose whether to save intermediate positions of particles
% - choose how many recursions can be done to simulate secondaries
% Note that choosing recursionLimit=0 means no secondaries will be modelled
simopts = SimulationOptions(exclude="all",threads=8,savesPositions=false,recursionLimit=2);

% choose which physical processes to model
% The list of supported models is given by SimulationOptions.allowedModels 
simopts.include("Ionisation","Bremsstrahlung","Nuclear","Bethe","MCS","Rutherford")
%% Run simulation
mcout = mcsimulate(mcinput,simopts);
%% Inspect results

% WARNING
% this produces LOTS of figures, please comment out uninteresting functions

% view the input geometry and generated samples
mcinput.view

% view the particle tracks inside the geometry
% - only consists of initial and final positions if savesPositions is false
mcout.viewTracks

% view the changes in transverse position and momentum for primaries
mcout.viewChanges

% view where the scattered particles reach the beam tube surface:
% - tube is representad as a 2D net with the y axis as its circumference
% - second plot will show volumetric energy deposition along same xy axes
% - "scalefactor" option scales the energy deposition by a factor
% - MaxThetaPoints is set to 20 to save computation time at higher samples
% - k means defines the # of regions for separate kernel density estimation
% Note that the KDE is inaccurate for most sample sizes under 1e5
mcout.viewIntersections("scalefactor",scaleFactor,"MaxThetaPoints",20,"kmeans",2);

% view the ranges of scattered particles in tungsten (W)
% - secondary ranges are of most interest, save time by excluding primaries
% - range estimation options are "CSDA", "True", or "Perpendicular"
mcout.viewCSDARange("W","method","Perpendicular")

