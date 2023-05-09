classdef Material < handle
    %MATERIAL A generic material with read-only properties
    %   This class obtains material properties from a formatted spreadsheet
    %   containing data for every material under consideration, and parses
    %   the spreadsheet to list its properties in a form MATLAB can
    %   interact with.
    %   The properties are all read only constants once the material object
    %   is initialised. Some of these properties will be vectors,
    %   corresponding to the properties of the individual elements that
    %   make up the material. Scalar properties represent the bulk behavior
    %   of the material. Units are listed next to the properties below.
    %   Properties that are not available will be saved as NaN and will
    %   throw a warning in the console
    %   To add additional properties relevent to your part of the project,
    %   extend this class, and draw additional data from the protected
    %   tabValues property (these will be as strings). You may then perform input
    %   validation on this data and save it as a property in the child
    %   class
    properties(GetAccess=protected,SetAccess=private)
        tabValues(1,:) table        % All material data from file stored as string
    end
    
    properties(SetAccess=protected)
        name(1,1) string            % Name of material, or chemical formula
        A(1,:) double               % AMU
        Z(1,:) double               % Dimensionless
        I(1,:) double               % MeV
        wt(1,:) double              % Dimensionless
        rho(1,1) double             % kg m^-3
        X0(1,1) double              % kg m^-2
        Tmax(1,1) double            % K
        isConductor(1,1) logical    % true if material has free conduction band electrons
    end

    properties(Dependent,SetAccess=protected)
        n(1,1) double               % number density of atoms (m^-3)
        x(1,1) double               % mole fractions of each element
        braggI(1,1) double          % approximate overall mean ionisation energy (MeV)
        
    end

    properties(Constant)
        % if a new materials spreadsheet is made, change this value
        defaultFilePath = "Material Properties Master.xlsx"

        % to use a different default material, change this
        defaultMaterial = "Chromox"
    end
    methods % Constructor
        function obj = Material(materialName,filename)
            %MATERIAL Construct an instance of this class
            %   If no filename is provided, the default name from Teams
            %   will be used
            if nargin == 0
                materialName = Material.defaultMaterial;
                filename = Material.defaultFilePath;
            elseif nargin == 1
                filename = Material.defaultFilePath;
            end
            matProps = Material.materialParser(filename,materialName);

            obj.A = Material.str2vec(matProps.AtomicWeight);
            obj.Z = Material.str2vec(matProps.AtomicNumber);
            obj.I = Material.str2vec(matProps.MeanIonisationEnergy);
            obj.wt = Material.str2vec(matProps.WeightPercentage);
            radLength = Material.str2vec(matProps.RadiationLength);
            if length(radLength)==1 % is bulk radiation length provided?
                obj.X0 = radLength;
            else
                obj.X0 = sum(obj.wt'.*radLength); % use approximation from Particle Physics Review 2012
            end
            obj.rho = Material.str2dub(matProps.Density);
            obj.Tmax = Material.str2dub(matProps.MaximumOperatingTemperature);
            obj.isConductor = logical(Material.str2dub(matProps.isConductor));
            obj.name = materialName;
        end
    end

    methods
        function E = dEdx(obj,particleBeta,particleMass,particleCharge)
            arguments
                obj(1,1) Material
                particleBeta(1,:) double = sqrt(1-1/(Consts.initialFourMomentum(1)/Consts.mpgev)^2) % v/c for incident particle
                particleMass(1,:) double = Consts.mpgev % rest mass of incident particle in GeV/c^2, defaults to proton mass
                particleCharge(1,:) double = 1;
            end
            K = Consts.Kbb; % GeV kmol^-1 m^2
            gamma = 1./sqrt(1-particleBeta.^2);
            me = Consts.megev*1000; % GeV/c^2
            M = particleMass*1000; % GeV/c^2
            Wmax = 2*me.*(gamma.^2-1)./(1+2.*gamma.*(me./M)+ (me./M).^2); % MeV
            logterm = 0.5.*log((2.*me.*gamma.^2.*particleBeta.^2.*Wmax)./(obj.braggI^2));
            ZoverA = sum(obj.wt.*obj.Z./obj.A)./Consts.MolarMassConstant; % kmol kg^-1

            % --------------delta term----------------
            IeV = obj.braggI*1e6; % eV
            rhogcm = obj.rho*1e-3; % g cm^-3
            plasmaEnergy = 28.816*(rhogcm*ZoverA)^.5; % eV
            xVal = log10(particleBeta.*gamma);
            C = 2*log(IeV/plasmaEnergy) + 1;
            if IeV < 100
                x1 = 2.0;
                if C < 3.681
                    x0 = .2;
                else
                    x0 = .326*C - 1.0;
                end
            else
                x1 = 3.0;
                if C < 5.215
                    x0 = .2;
                else
                    x0 = .326*C - 1.5;
                end
            end
            delta = zeros(1,length(gamma));
            delta(xVal>=x1) = 2*log(10).*xVal(xVal>=x1) - C;
            delta(xVal>=x0) = 2*log(10)*xVal(xVal>=x0) - C + ((C - 2*log(10)*x0)/((x1-x0)^3)).*(x1-xVal(xVal>=x0)).^3;
            % remaining delta vals are zero
            % -----------------------------------------------
            Iproton = (particleCharge==1 & particleMass==Consts.mpgev);
            E = zeros(1,length(particleMass));
            E(Iproton) = obj.rho*K.*particleCharge(Iproton).^2*ZoverA./particleBeta(Iproton).^2.*(logterm(Iproton) -particleBeta(Iproton).^2 -delta(Iproton)./2); % MeV
            
            E(~Iproton) = obj.rho*K*ZoverA./particleBeta(~Iproton).^2 .* (0.5*log((2.*(gamma(~Iproton)+1))/(obj.braggI/(1e3*Consts.megev))^2) + 0.5.*F(particleBeta(~Iproton),gamma(~Iproton)-1,0.5*(gamma(~Iproton)-1))-delta(~Iproton)./2);
            % add in effects of brehmstrahlung for electrons
            E(~Iproton) = E(~Iproton) + (gamma(~Iproton)-1).*particleMass(~Iproton)./(obj.X0./obj.rho);
        end
    
        function delX = csdaRange(obj,particleKineticEnergy,particleMass,particleCharge)
            % determine limit of validity
            mke = Consts.MinimumKineticEnergy(particleMass,particleCharge);

            % define function for ode
            function val = dxdE(kineticEnergy,~)
                    totalEnergy = kineticEnergy + particleMass;
                    gamma = totalEnergy./particleMass;
                    beta = sqrt(1-1./gamma.^2);
                    val = -1./obj.dEdx(beta,repmat(particleMass,1,length(beta)),repmat(particleCharge,1,length(beta)));
            end

            % perform csda calculation
            if particleKineticEnergy<=mke
                % just use the low energy range
                delX = obj.lowEnergyCSDA(particleKineticEnergy,particleMass,particleCharge);
            else % integrate and then apply low energy
                [~,delX] = ode78(@dxdE,[particleKineticEnergy mke],0); % integrate to mke
                delX = delX(end) + obj.lowEnergyCSDA(mke,particleMass,particleCharge); % add low energy correction
            end
        end

        function Efinal = energyChange(obj,particleKineticEnergy,particleMass,particleCharge,deltaX)
            arguments
                obj(1,1) Material
                particleKineticEnergy(:,1) double {mustBeNonnegative}
                particleMass(:,1) double {mustBeNonnegative}
                particleCharge(:,1) double {mustBeInteger}
                deltaX(:,1) double {mustBePositive}
            end
            function val = dEdx(~,kineticEnergy)
                    totalEnergy = kineticEnergy + particleMass;
                    gamma = totalEnergy./particleMass;
                    beta = sqrt(1-1./gamma.^2);
                    val = -1*obj.dEdx(beta,particleMass,particleCharge);
                    val = val';
            end
            
            [~,E] = ode23(@dEdx,[0 min(deltaX)],particleKineticEnergy);
            Efinal = E(end,:);
        end
    end

    methods % getters
        function n = get.n(obj)
            n = Consts.Na*obj.rho*sum(obj.wt./obj.A);
        end

        function I= get.braggI(obj)
            numerator = sum(obj.wt.*(obj.Z./obj.A).*log(obj.I));
            denominator = sum(obj.wt.*(obj.Z./obj.A));
            I = exp(numerator/denominator);
        end

        function x = get.x(obj)
            totalMoles = sum(obj.wt./obj.A);
            x = obj.wt./(obj.A.*totalMoles);
        end


    end

    methods(Access=protected)
        function r = lowEnergyCSDA(obj,particleKineticEnergy,particleMass,particleCharge)
            if particleMass == Consts.mpgev && particleCharge == 1
                % particle is proton
                % use Bortfeld relation
                alpha = 2.2e-3; % cm per MeV
                p = 1.77;
                % requires energy in MeV and gives result in cm
                rcm = alpha * (particleKineticEnergy*1000)^p;
                % convert to metres
                r = rcm/100;

                % this is the range in water, so needs scaling
                % proportional to effective atomic mass of compound
                % inv proportional to density
                AeffWater = 14.4444;
                rhoWater = 1000;
                r = sum(obj.wt.*obj.A)./AeffWater * rhoWater./obj.rho * r;
            elseif particleMass == Consts.megev && particleCharge == -1
                % particle is electron
                % Use D. Tan relation
                % fit is for r in g/cm^2 and Ek in keV
                rgcm = 1.90e-6*(particleKineticEnergy*1e6)^1.6 * sum((obj.A./obj.Z).^2.5 .*obj.wt);
                % convert to kg m^-2, then divide by density
                r = rgcm*10./obj.rho;
            else
                r=0;
            end
        end
    end

    methods(Static,Access=protected) % Parser functions
        function materialProps = materialParser(filename,material)
            %MATERIALPARSER Parses the provided material spreadsheet for a
            % specific material
            %   The table is read as a string so that multiple
            %   comma-separated values can be held in each spreadsheet
            %   cell. This also allows for the inclusion of non-numeric
            %   data types in the future if needed
            %   Inputs:
            %       filename        -The file location of the spreadsheet
            %                       containing the material data. It should
            %                       be saved in Excel (xlsx,xls) format.
            %       material        -The name of the material
            %                       (case-sensitive) for the data to be
            %                       drawn from
            %   Outputs:
            %       materialProps   -Table with one row corresponding to
            %                       the properties of the provided
            %                       material. All properties are stored as
            %                       strings, so later input validation is
            %                       needed to convert these to correct data
            %                       types.
            opts = detectImportOptions(filename);
            opts = opts.setvartype("string");
            materialProps = readtable(filename,opts);
            materialProps = materialProps(materialProps.Material==material,:);
            if isempty(materialProps)
                throw(MException("Material:InvalidMaterial","The material name provided could not be found in the spreadsheet. Note that the parser is case-sensitive"))
            end
        end
        function v = str2vec(s,NegativeNaN)
            arguments
                s(1,1) string
                NegativeNaN(1,1) logical = true
            end
            %STR2VEC converts a comma-separated string into a vector of
            %doubles
            %   Inputs:
            %       s               -String with comma-separated numerical
            %                       values
            %   Outputs:
            %       v               -Double vector from the input string,
            %                       NaN if string is invalid
            try
                vc = textscan(s,'%f','Delimiter',',',"ReturnOnError",false);
                v = vc{1};
                if (any(v<0)&&NegativeNaN)
                    v = NaN;
                end
            catch me
                v = NaN;
            end
            if any(isnan(v))
                warning("Some properties are NaN")
            end
        end
        function x = str2dub(s,NegativeNaN)
            arguments
                s(1,1) string
                NegativeNaN(1,1) logical = true
            end
            try
                x = double(s);
                if (x<0&&NegativeNaN)
                    x = NaN;
                end
            catch me
                x = NaN;
            end
            if isnan(x)
                warning("Some properties are NaN")
            end
        end
    end

    methods(Static) % static, spreadsheet overviews
        function mats = allNames(filename)
            arguments
                filename(1,1) string = Material.defaultFilePath
            end
            opts = detectImportOptions(filename);
            opts = opts.setvartype("string");
            materialProps = readtable(filename,opts);
            mats = sort(materialProps.Material);
        end
    end
end

function f = F(beta,tau,tauUp) % Berger-Seltzer correction (GEANT4)
    gamma = 1./sqrt(1-beta.^2);
    f = -1 -beta.^2 +log((tau-tauUp).*tauUp) + tau./(tau-tauUp) + 1./gamma.^2 .* (tauUp.^2./2 + (2.*tau+1).*log(1-tauUp./tau));
end