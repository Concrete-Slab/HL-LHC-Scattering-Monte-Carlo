function val = runSim(material,simopts,thickness,sideLength,offset,pdf,expectedProtons,nSamples)
%     log10nMax = 7;
%     log10nMin = 4;
%     log10tMax = 6;
%     log10tMin = 2;
%     nSamples = ceil(10^((log10nMax-log10nMin)./(log10tMax-log10tMin) * (-log10(thickness)-log10tMin)+log10nMin));
%     nSamples = ceil(10^(-2/3*log10(thickness)+8/3));
    scaleFactor = expectedProtons/nSamples;
    geometry = OffsetRectangle(sideLength,sideLength,offset,thickness,0);
    geometry.pdf = pdf;
    samples = geometry.generate(nSamples,"echo","off");
    in = MCInput(geometry,material,samples=samples);
    out = mcsimulate(in,simopts,"echo","off");
    [~,~,E] = out.energyDistribution("MaxThetaPoints",20,"echo","off");
    E = E.*scaleFactor.*1e-3;
    val = max(max(E));
end