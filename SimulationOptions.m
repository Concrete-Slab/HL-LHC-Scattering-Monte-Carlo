classdef SimulationOptions < handle
    %SIMULATIONOPTIONS Summary of this class goes here
    %   To add a process, add its string name to the softModels or
    %   hardModels constant array (or add to hidden ones), and then add its
    %   constructor to the switch-case in getSoft() or getHard()
    %   respectively. If the process produces secondaries (and is hard),
    %   also add it to the produceSecondaries hidden list

    properties(Constant)
%         allowedModels(1,4) string = ["MCS","Rutherford","Photoionisation","Bethe"]
        softModels(1,:) string = ["MCS","Bethe"]
        hardModels(1,:) string = ["Rutherford","Ionisation"]
        allowedModels(1,:) string = [SimulationOptions.softModels,SimulationOptions.hardModels]
    end

    properties(Hidden,Constant)
%         allModels = ["MCS","Rutherford","Photoionisation","Bethe","DummyHardProcess","DummySoftProcess"]
        hiddenModels(1,:) string = ["DummyHardProcess","DummySoftProcess","DummySecondaryProcess"]
        allSoft(1,:) string = [SimulationOptions.softModels,"DummySoftProcess"]
        allHard(1,:) string = [SimulationOptions.hardModels,"DummyHardProcess","DummySecondaryProcess"]
        allModels(1,:) string = [SimulationOptions.allSoft,SimulationOptions.allHard];
        produceSecondaries(1,:) string = ["DummySecondaryProcess","Ionisation"]
        defaultMaximumSecondaries(1,1) double {mustBeNonnegative,mustBeInteger} = 100;
        defaultPoissonPercentile(1,1) double {mustBeLessThan(defaultPoissonPercentile,1),mustBePositive} = 0.9995
    end

    properties
        threads(1,1) double {mustBeInteger,mustBePositive,mustBeLessThanOrEqual(threads,8)} = 1
        savesPositions(1,1) logical = true
        recursionLimit(1,1) double {mustBeInteger,mustBeNonnegative} = 0;
        maximumSecondaries(1,1) double {mustBeNonnegative,mustBeInteger} = SimulationOptions.defaultMaximumSecondaries
        poissonPercentile(1,1) double {mustBeNonnegative,mustBeLessThan(poissonPercentile,1)} = SimulationOptions.defaultPoissonPercentile
    end

    properties(SetAccess=private,Hidden)
        excluded(1,:) string {SimulationOptions.mustBeModel}
    end

    properties(Dependent,SetAccess=private)
        included(1,:) string
    end

    methods % constructor
        function obj = SimulationOptions(options)
            %SIMULATIONOPTIONS Construct an instance of this class
            %   Detailed explanation goes here
            arguments
                options.exclude(1,:) string {SimulationOptions.mustBeModel} = [];
                options.threads(1,1) double {mustBeInteger,mustBePositive,mustBeLessThanOrEqual(options.threads,8)} = 1
                options.savesPositions(1,1) logical = true
                options.recursionLimit(1,1) double {mustBeInteger,mustBeNonnegative} = 0;
                options.poissonPercentile(1,1) double {mustBeNonnegative,mustBeLessThan(options.poissonPercentile,1)} = SimulationOptions.defaultPoissonPercentile;
                options.maximumSecondaries(1,1) double {mustBeNonnegative,mustBeInteger} = SimulationOptions.defaultMaximumSecondaries;
            end
            obj.exclude(options.exclude);
            obj.threads = options.threads;
            obj.savesPositions = options.savesPositions;
            obj.recursionLimit = options.recursionLimit;
            obj.poissonPercentile = options.poissonPercentile;
            obj.maximumSecondaries = options.maximumSecondaries;
        end
    end


    methods % getters
        function inc = get.included(obj)
            idx = ~ismember(SimulationOptions.allModels,obj.excluded);
            inc = SimulationOptions.allModels(idx);
        end
    end

    methods % public methods
        function include(obj,modelVector,model)
            arguments
                obj(1,1) SimulationOptions
                modelVector(1,:) string {SimulationOptions.mustBeModel}
            end
            arguments(Repeating)
                model(1,1) string {SimulationOptions.mustBeModel}
            end
            models = [modelVector string(model)];
            if ismember('all',models)
                obj.excluded = SimulationOptions.hiddenModels;
            else
                if ismember('hard',models)
                    models = unique([models,SimulationOptions.hardModels]);
                elseif ismember('soft',models)
                    models = unique([models,SimulationOptions.softModels]);
                end
                obj.excluded = obj.excluded(~ismember(obj.excluded,models));
            end
        end

        function exclude(obj,modelVector,model)
            arguments
                obj(1,1) SimulationOptions
                modelVector(1,:) string {SimulationOptions.mustBeModel}
            end
            arguments(Repeating)
                model(1,1) string {SimulationOptions.mustBeModel}
            end
            models = [modelVector string(model)];
            if ismember("all",models)
                models = SimulationOptions.allModels;
            elseif ismember("soft",models)
                models = unique([models,SimulationOptions.allSoft]);
            elseif ismember("hard",models)
                models = unique([models,SimulationOptions.allHard]);
            end
            models = models(~ismember(models,["soft","hard"]));
            nonmembers = ~ismember(models,obj.excluded);
            obj.excluded = [obj.excluded models(nonmembers)];
        end

        function [hardsim,softsim] = configure(obj,mcinput)
            % configure hard and soft simulation objects
            arguments(Input)
                obj(1,1) SimulationOptions
                mcinput(1,1) MCInput        % geometry and material info.
            end
            arguments(Output)
                hardsim(1,:) HardProcess
                softsim(1,:) SoftProcess
            end
            %             if isempty(obj.included)
            %                 error("No included models!")
            %             end
            softStr = unique(SimulationOptions.allSoft(ismember(SimulationOptions.allSoft,obj.included)));
            hardStr = unique(SimulationOptions.allHard(ismember(SimulationOptions.allHard,obj.included)));
            
            if ismember("Ionisation",hardStr)
                % need to configure a bethe bloch process alongside
                if ismember("Bethe",softStr)
                    nSoft = length(softStr);
                else
                    nSoft = length(softStr) + 1;
                end
                hardsim(length(hardStr)) = Ionisation(mcinput.material);
                softsim(nSoft) = SoftBetheBloch(mcinput.material,Consts.data.ionisationTcut);
                hardStr = hardStr(~ismember(hardStr,"Ionisation"));
                softStr = softStr(~ismember(softStr,"Bethe"));
            else
                hardsim = HardProcess.empty(length(hardStr),0);
                softsim = SoftProcess.empty(length(softStr),0);
            end
            % TODO preallocate
            for isoft = length(softStr):-1:1
                softsim(isoft) = getSoft(softStr(isoft),mcinput);
            end
            for ihard = length(softStr):-1:1
                hardsim(ihard) = getHard(hardStr(ihard),mcinput);
            end
        end

        function n = estimatedSecondaries(obj,mcinput)
            arguments(Input)
                obj(1,1) SimulationOptions
                mcinput(1,1) MCInput
            end
            arguments(Output)
                n(1,1) double {mustBeNonnegative,mustBeInteger}
            end
            idx = and(ismember(obj.included,SimulationOptions.allHard),ismember(obj.included,SimulationOptions.produceSecondaries));
            n=0;
            defaultEnergy = Consts.initialFourMomentum(1);
            defaultCharge = 1;
            defaultMass = Consts.mpgev;
            if obj.recursionLimit>0
                lambdaTot = 0;
                hardstr = obj.included(idx);
                for i=1:length(hardstr)
                    hard = getHard(hardstr(i),mcinput);
                    hard.particle = ParticleHandle(defaultEnergy,defaultCharge,defaultMass);
                    lambdaTot = lambdaTot + hard.lambda;
                    hard.particle.delete
                end
                meanInteractions = mcinput.geometry.thickness/lambdaTot;
                % # interactions ~ poisson distribution with above mean
                % inverse poisson cdf solver
                P = 0;
                n = 0;
                while P<obj.poissonPercentile
                    P = P + exp(-meanInteractions)*meanInteractions^n./gamma(n + 1);
                    n = n + 1;
                end
                if n>obj.maximumSecondaries
                    warning("Predicted number of secondaries is %d, which is above the maximum limit of %d, accuracy may be reduced",n,obj.maximumSecondaries)
                end
                n = min(n,obj.maximumSecondaries);
            end
        end
    end

    methods(Static,Access=private) % input validation
        function mustBeModel(model)
            mustBeMember(model,[SimulationOptions.allModels,'all','hard','soft'])
        end
    end
end

function sp = getSoft(strp,inp)
arguments(Input)
    strp(1,1) string
    inp(1,1) MCInput
end
arguments(Output)
    sp(1,:) SoftProcess {mustBeScalarOrEmpty}
end
switch strp
    case "MCS"
        sp = MCS(inp.material);
    case "DummySoftProcess"
        sp = DummySoftProcess;
    case "Bethe"
        sp = SoftBetheBloch(inp.material);
    otherwise
        warning("Couldnt find soft model %s, returning empty",strp)
        sp = SoftProcess.empty(1,0);
end
end

function hp = getHard(strp,inp)
arguments(Input)
    strp(1,1) string
    inp(1,1) MCInput
end
arguments(Output)
    hp(1,:) HardProcess {mustBeScalarOrEmpty}
end
switch strp
    case "DummyHardProcess"
        hp = DummyHardProcess(inp.material);
    case "DummySecondaryProcess"
        hp = DummySecondaryProcess(inp.material);
    case "Ionisation"
        hp = Ionisation(inp.material);
    case "Rutherford"
        hp = RutherfordScattering(inp.material,inp.geometry);
    otherwise
        warning("Couldnt find hard model %s, returning empty",strp)
        hp = HardProcess.empty(1,0);
end
end