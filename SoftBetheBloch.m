classdef SoftBetheBloch < SoftProcess
    %SOFTBETHEBLOCH Summary of this class goes here
    %   Detailed explanation goes here

    properties(SetAccess=private)
        I(1,1) double               % bragg mean ionisation energy (MeV)
        Kextended(1,1) double       % material constant GeV m^-1
        x1(1,1) double              % delta correction term
        x0(1,1) double              % delta correction term
        C(1,1) double               % delta correction term
        Tcut(1,1) double            % minimum energy transfer MeV
        nel(1,1) double             % number density of electrons, m^-3
        isConductor(1,1) logical    % is material a conductor
    end

    properties(SetAccess=private,Dependent)
        delta(1,1) double           % delta correction
        Tmax(1,1) double            % maximum energy transfer (MeV)
        x(1,1) double               % x value for delta correction
        delta0(1,1) double          % delta correction term
        Tup(1,1) double             % maximum of Tcut and Tmax
    end

    properties(Dependent,SetAccess=protected)
        dEdx % GeV m^-1
    end

    properties(Constant)
        K(1,1) double = Consts.Kbb % GeV kmol^-1 m^2
    end

    methods % constructor
        function obj = SoftBetheBloch(material,Tcut,particle)
            %SOFTBETHEBLOCH Construct an instance of this class
            %   Detailed explanation goes here
            arguments
                material(1,1) Material = Material
                Tcut(1,1) double = Inf      % cutoff energy using hard ionisation, MeV
                particle(1,:) ParticleHandle {mustBeScalarOrEmpty} = ParticleHandle.empty
            end
            obj@SoftProcess(particle)
            % FIND VALUES FOR DELTA0
            obj.isConductor = material.isConductor;
            ZoverA = sum(material.wt.*material.Z./material.A)./Consts.MolarMassConstant; % kg kmol^-1
            rho = material.rho; % kg m^-3
            obj.I = material.braggI; % MeV
            obj.Kextended = Consts.Kbb.*ZoverA*rho; % GeV m^-1
            plasmaEnergy = 28.816*(rho*1e-3*ZoverA)^.5; % eV
            obj.C = 2*log(obj.I*1e6/plasmaEnergy) + 1; % ln(eV/eV)
            if obj.I*1e6 < 100
                obj.x1 = 2.0;
                if obj.C < 3.681
                    obj.x0 = .2;
                else
                    obj.x0 = .326*obj.C - 1.0;
                end
            else
                obj.x1 = 3.0;
                if obj.C < 5.215
                    obj.x0 = .2;
                else
                    obj.x0 = .326*obj.C - 1.5;
                end
            end
            obj.Tcut = Tcut; % MeV
            obj.nel = Consts.Na*sum(material.wt.*material.Z./material.A)*material.rho; % m^-3
        end

    end

    methods % getters
        function e = get.dEdx(obj)
            arguments(Output)
                e(1,1) double % dEdx, GeV m^-1
            end
            gamma = obj.particle.gamma;
            beta = sqrt(1-1/gamma^2);
            if obj.particle.mass == Consts.megev && obj.particle.charge == -1
                % electron, use slightly different correlation
                % Berger-Seltzer formula
                tau = gamma-1;
                me = Consts.megev*1e3;
                tauUp = min(tau/2,obj.Tcut./me);
                fterm = F(tau,tauUp);
                logterm = log((2*(gamma+1))/(obj.I/me)^2);
                e = -2*pi*Consts.re^2*Consts.megev*obj.nel/beta^2 * (logterm + fterm - obj.delta);
            elseif obj.particle.mass == Consts.mpgev && obj.particle.charge == 1
                % proton, use corrected bethe bloch formula for T<Tcut
                gamma = obj.particle.gamma;
                z = obj.particle.charge;
                tmax = obj.Tmax;
                tup = min(obj.Tcut,tmax);
                % spin correction
                S = (0.5*tup/(obj.particle.energy*1e3))^2;
                e = -obj.Kextended*z^2/beta^2 * (0.5*log(2*Consts.megev*1e3*(gamma^2*beta^2)*obj.Tmax / obj.I^2) - 0.5*beta^2*(1+tup/tmax) - obj.delta/2 + S/2);
            else
                e = 0;
            end
        end
        function tmax = get.Tmax(obj)
            if obj.particle.charge ==-1 && obj.particle.mass==Consts.megev
                % electron
                tmax = 0.5*Consts.megev*1e3*(obj.particle.gamma-1); % MeV
            else
                % hadron
                gamma = obj.particle.gamma;
                me = Consts.megev;
                M = obj.particle.mass;
                tmax = 2*me*1e3*(gamma^2-1)./(1+2*gamma*(me/M)+(me/M)^2); % MeV
            end
        end

        function del = get.delta(obj)
            xval = obj.x;
            if xval>=obj.x1
                del = 2*log(10)*xval-obj.C;
            elseif xval>=obj.x0
                t1 = obj.C - 2*log(10)*obj.x0;
                t2 = (obj.x1-obj.x0)^3;
                a = t1/t2;
                del = 2*log(10)*xval - obj.C + a*(obj.x1-xval);
            elseif xval<obj.x0 && obj.isConductor
                del = obj.delta0*10^(2*xval-obj.x0);
            else
                del = 0;
            end
        end

        function xVal = get.x(obj)
            xVal = log10(obj.particle.gamma^2-1);
        end

        function d0 = get.delta0(obj)
            d0 = 4.606*obj.x-obj.C;
        end

        function tup = get.Tup(obj)
            tup = min(obj.Tcut*1e3,obj.Tmax);
        end
    end

    methods % inherited
        function [dE,dAngle,dPos] = update(obj,delX)
            arguments(Input)
                obj(1,1) SoftBetheBloch
                delX(1,1) double
            end
            arguments(Output)
                dE(1,1) double
                dAngle(2,1) double
                dPos(2,1) double
            end
            % bethe bloch equation integrated at constant value for delX
            dE = -delX * obj.dEdx;
            dAngle = [0;0];
            dPos = [0;0];
        end
    end
end

function f = F(tau,tauUp)
    % GEANT4 handbook page 122 (136)
    gamma = tau + 1;
    beta = sqrt(1-1/gamma^2);
    f = -1 -beta^2 + log((tau-tauUp)*tauUp) + tau/(tau-tauUp) + 1/gamma^2 * (tauUp^2/2 + (2*tau+1)*log(1-tauUp/tau));
end
