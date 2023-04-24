classdef RuntimeData < handle
    %RUNTIMEDATA Summary of this class goes here
    %   After adding a new variable, specify a default value under the constant
    %   properties header, and add a line setting your variable to default
    %   under the default(obj) method.
    
    properties
        % Tcut for ionisation model, MeV
        ionisationTcut(1,1) double {mustBePositive} = RuntimeData.defaultIonisationTcut
        % gaussian fit model for multiple coulomb scattering
        mcsFit(1,1) string {mustBeMember(mcsFit,["LynchDahl","Highland"])} = RuntimeData.defaultMCSFit
        % number of points used to integrate the pdf over a geometry
        cdfResolution(1,1) double {mustBeInteger,mustBePositive} = RuntimeData.defaultCDFResolution
        % suggested value for the grid density parameter of nonUniformGrid
        gridDensityMultiplier(1,1) double {mustBePositive} = RuntimeData.defaultgridDensityMultiplier
        % suggested function to use for non-uniform grid generation
        gridDensityFunction(1,1) function_handle = RuntimeData.defaultgridDensityFunction
        % bunch spacing, seconds
        bunchSpacing(1,1) double {mustBeMember(bunchSpacing,[25e-9,50e-9])} = RuntimeData.defaultBunchSpacing
    end

    properties(Constant,Hidden)
        defaultIonisationTcut = 2; % MeV
        defaultMCSFit = "LynchDahl";
        defaultCDFResolution = 100;
        defaultgridDensityMultiplier = 10;
        defaultgridDensityFunction(1,1) function_handle = @sinh;
        defaultBunchSpacing = 25e-9;
    end
    
    methods % public
        function obj = RuntimeData
            %RUNTIMEDATA Construct an instance of this class
            %   Detailed explanation goes here
            default(obj)
        end

        function default(obj)
            obj.ionisationTcut = obj.defaultIonisationTcut;
            obj.mcsFit = obj.defaultMCSFit;
            obj.cdfResolution = obj.defaultCDFResolution;
            obj.gridDensityMultiplier = obj.defaultgridDensityMultiplier;
            obj.gridDensityFunction = obj.defaultgridDensityFunction;
            obj.bunchSpacing = obj.defaultBunchSpacing;
        end
    end

    properties(Dependent,SetAccess=private) 
        Nb % protons per bunch
    end

    methods % dependent property getters
        function n = get.Nb(obj)
            if obj.bunchSpacing==25e-9
                n = 2.2e11;
            elseif obj.bunchSpacing==50e-9
                n = 3.5e11;
            else
                n = 1.15e11;
            end
        end
    end

end

% add dcf grid density + function

