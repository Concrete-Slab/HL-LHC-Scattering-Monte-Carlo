classdef Ionisation < HardProcess
    %IONISATION Summary of this class goes here
    %   Detailed explanation goes here


    properties(Dependent)
        sigmaTot            % m^2
        Tmax(1,1) double    % maximum energy transfer in GeV
    end

    properties(SetAccess=private)
        Z(1,1) double       % mean atomic number
        Tcut(1,1) double    % minimum energy transfer in GeV
        I(1,1) double       % overall mean ionisation energy (bragg) GeV
    end


    methods
        function obj = Ionisation(material)
            %IONISATION Scattering of primary particles by ionisation of a
            %delta ray. Subclass of HardProcess.
            %   GEANT4 Physics documentation page 175 gives the governing
            %   equations for the behaviour of this scattering process
            arguments
                material(1,1) Material
            end
            obj@HardProcess(material);
            obj.I = material.braggI/1e3;
            obj.Z = sum(material.x.*material.Z); % average by mole fraction rather than weight
            % find a better description of Tcut, not just = I
            obj.Tcut = Consts.data.ionisationTcut/1e3; % GeV
        end

        function tmax = get.Tmax(obj)
            me = Consts.megev;
            gamma = obj.particle.gamma;
            if obj.particle.mass==me && obj.particle.charge == -1
                % electron transfers at most half of its kinetic energy
                tmax = 0.5*me*(gamma-1);
            else
                % formula for other particles...
                tmax = 2*me*(gamma^2-1)./(1+2*gamma*(me/obj.particle.mass)+ (me/obj.particle.mass)^2);
            end
        end

        function st = get.sigmaTot(obj)
            % four momentum in frame of reference where primary travels
            % along the z axis
            beta = obj.particle.beta;
            E = obj.particle.energy; % GeV
            z = obj.particle.charge; % e
            spin = 1;
            me = Consts.me*Consts.c^2/(1e9*Consts.e); % GeV/c^2
            tmax = obj.Tmax; % GeV
            % cross section in m^2
            if obj.Tcut>tmax
                st=0;
            elseif obj.particle.charge == 1 && obj.particle.mass==Consts.mpgev
                % particle is a proton
                st = 2*pi*Consts.re^2*obj.Z*z^2*me/beta^2 * ((1/obj.Tcut-1/tmax)-beta^2/tmax*log(tmax/obj.Tcut)+spin*(tmax-obj.Tcut)/(2*E^2));
            elseif obj.particle.charge == -1 && obj.particle.mass==Consts.megev
                % particle is an electron
                x = obj.Tcut./(obj.particle.energy-obj.particle.mass);
                gamma = obj.particle.gamma;
                beta = sqrt(1-1/gamma^2);
                st = 2*pi*Consts.re^2*obj.Z/(beta^2*(gamma-1)) * ((gamma-1)^2/gamma^2 * (0.5-x) + 1/x - 1/(1-x) - (2*gamma-1)/gamma^2 * log((1-x)/x));
            else
                st = 0;
            end
        end

        function [dE,dAngles,secondaryParticle] = interact(obj)
            arguments(Input)
                obj(1,1)
            end
            arguments(Output)
                dE(1,1) double
                dAngles(4,1) double
                secondaryParticle(6,:) double
            end
%             acceptedValue = false;
%             tmax = obj.Tmax;
%             T=0;
%             % TRY VARIABLE PRECISION ARITHMETIC?

%             while ~acceptedValue && tmax>obj.Tcut
%                 z = rand(2,1);
%                 % sample T from normalised f(t) = 1/T^2
%                 T = 1/(1/obj.Tcut - z(1)*(1/obj.Tcut-1/tmax));
%                 % obtain g(T)
%                 mingT = obj.gProton(tmax);
%                 maxgT = obj.gProton(obj.Tcut);
%                 gTnormalised = (obj.gProton(T)-mingT)/(maxgT-mingT);
%                 if z(2)<gTnormalised
%                     acceptedValue=true;
%                     % Some of the energy transfer goes towards exciting the atom
%                     T = sym(T-obj.I);
%                 end
%             end
            E = sym(obj.particle.energy); % GeV symbolic
            m = sym(obj.particle.mass); % GeV/c^2
            if obj.particle.mass == Consts.mpgev && obj.particle.charge == 1
                % particle is proton
                T = sym(obj.protonEvaluate);
            elseif obj.particle.mass == Consts.megev && obj.particle.charge == -1
                % particle is electron
                T = sym(obj.electronEvaluate);
            else
                T = sym(0);
            end
            % we now know the energy transfer to the secondary particle.
            % need to use energy-momentum conservation to convert this to
            % convert this to a polar angle theta. azimuth is isotropic
            me = sym(Consts.megev); % GeV
            % energy conservation
            p1p = sqrt(E^2-m^2);
            p2p = sqrt((E-T)^2-m^2);
            p2e = sqrt((me+T)^2-me^2);
            % momentum conservation
            theta = acos((p1p^2+p2p^2-p2e^2)/(2*p1p*p2p));


            cosTheta = cos(theta);
            sinTheta = sin(theta);
            dE = -double(T);
%             obj.particle.energy = double(E-T);
            % azimuthal angle (phi) is asotropic
            phi = 2*pi*rand;
            dAngles=double([sin(phi);cos(phi);sinTheta;cosTheta]);
            if any(~isreal(dAngles))
                dAngles
            end
            if dE<0
                secondaryMass = double(me);
                secondaryCharge = -1;
                secondaryEnergy = secondaryMass+double(T);
%                 secondaryMomentumMagnitude = p2e;
%                 sinSecondaryTheta = -sin(theta)*p2p/p2e;
%                 cosSecondaryTheta = (initialMomentum-finalMomentum*cosTheta)/secondaryMomentumMagnitude;
                secondaryPhi = phi+pi;
%                 secondaryMomentum = secondaryMomentumMagnitude.*[cos(secondaryPhi)*sinSecondaryTheta;sin(secondaryPhi)*sinSecondaryTheta;cosSecondaryTheta];
                secondaryMomentum = double([cos(secondaryPhi)*(-sin(theta)*p2p);sin(secondaryPhi)*(-sin(theta)*p2p);p2e*cos(theta)]);
                secondaryParticle = [secondaryEnergy;double(secondaryMomentum);secondaryCharge;secondaryMass];
            else
                secondaryParticle = double.empty(6,0);
            end
        end
    end

    methods(Access=private)
        function gT = gProton(obj,T)
            beta = obj.particle.beta;
            tmax = obj.Tmax;
            E = obj.particle.energy;
            gT = 1-beta^2*T/tmax + T^2/(2*E^2);
        end

        function ge = gElectron(obj,e)
            arguments(Output)
                ge(1,1) double
            end
            gamma = obj.particle.gamma;
            ge = 4/(9*gamma^2-10*gamma+5)*((gamma-1)^2*e^2 - (2*gamma^2 + 2*gamma -1)*e/(1-e) + gamma^2/((1-e)^2));
        end

        function T = electronEvaluate(obj)
            acceptedValue = false;
            e0 = obj.Tcut/(obj.particle.energy-obj.particle.mass);
            gmax = obj.gElectron(e0);
            gmin = obj.gElectron(0.5);
            T=0;
            while ~acceptedValue
                z = rand(2,1);
                e = 1/(1/e0-(1/e0-2)*z(1));
                ge = (obj.gElectron(e)-gmin)./(gmax-gmin);
                if z(2) < ge
                    acceptedValue = true;
                    T = e*(obj.particle.energy-obj.particle.mass)-obj.I;
                end
            end
            
        end

        function T = protonEvaluate(obj)
            acceptedValue = false;
            tmax = obj.Tmax;
            T = 0;
            mingT = obj.gProton(tmax);
            maxgT = obj.gProton(obj.Tcut);
            while ~acceptedValue && tmax>obj.Tcut
                z = rand(2,1);
                % sample T from normalised f(t) = 1/T^2
                T = 1/(1/obj.Tcut - z(1)*(1/obj.Tcut-1/tmax));
                % obtain g(T)
                gTnormalised = (obj.gProton(T)-mingT)/(maxgT-mingT);
                if z(2)<gTnormalised
                    acceptedValue=true;
                    % Some of the energy transfer goes towards exciting the atom
                    T = T-obj.I;
                end
            end
        end
    end
end

