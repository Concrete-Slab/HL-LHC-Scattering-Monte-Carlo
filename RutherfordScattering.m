classdef RutherfordScattering <HardProcess
    %RUTHERFORDSCATTERING Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        K(1,1) double           % constant = 4*pi*alpha^2*(hbarc)^2*Z^2
        MN(1,1) double          % average nuclear mass GeV/c^2
        Rsquared(1,1) double    % average nuclear radius squared, m^2
        maxDeltaX(1,1) double   % thickness of material or the maximum step size
        Zi(1,:) double          % atomic numbers of material
        Ai(1,:) double          % atomic masses of material
        X0(1,1) double          % radiation length kg m^-2
        rho(1,1) double         % density of material kg m^-3
        wt(1,:) double          % weight percentages of each element
    end

    properties(Dependent)
        sigmaTot                % m^2
        tmin(1,1) double        % GeV^2
        
    end
    
    methods
        function obj = RutherfordScattering(material,geometry,particle)
            %RUTHERFORDSCATTERING Construct an instance of this class
            %   Detailed explanation goes here
            arguments
                material(1,1) Material = Material
                geometry(1,1) Geometry = OffsetRectangle(0.01,0.01,0.01,0.001,0);
                particle(1,:) ParticleHandle {mustBeScalarOrEmpty} = ParticleHandle.empty;
            end
            obj@HardProcess(material,particle);
            obj.maxDeltaX = geometry.thickness/(cos(geometry.xRotation*pi/180));
            A = sum(material.wt.*material.A);
            obj.Rsquared = (1.2e-15 * A^0.33)^2;
            obj.MN = (A*Consts.MolarMassConstant./Consts.Na) * Consts.c^2/(1e9*Consts.e);
            obj.K = 4*pi*Consts.alpha^2*(Consts.hbarc/(1e9*Consts.e))^2*sum(material.wt.*material.Z.^2);
            obj.Zi = material.Z;
            obj.Ai = material.A;
            obj.X0 = material.X0;
            obj.rho = material.rho;
            obj.wt = material.wt;
        end

        function st = get.sigmaTot(obj)
            % TODO Re-evaluate this expression
            if obj.particle.charge~=1
                st = 0;
            else
                % integrate K*1/t^2*exp(-0.856e3 *t*R^2) from tmin to
                % infinity
                % integral evaluates to statement below
                tm = obj.tmin;
                a = 0.856e3*obj.Rsquared;
                % integration by parts, second part evaluates to upper
                % incomplete gamma function
                st = obj.K*((exp(-a*tm)/tm)-gammainc(0,a*tm,'upper'));
            end
        end

        function t = get.tmin(obj)
%             t = 0.998e-3;
            % change this so it corresponds to minimum scattering angle
            theta0 = MCS.getTheta0(obj.particle.momentum,obj.particle.beta,obj.particle.charge,obj.Zi,obj.Ai,obj.X0,obj.rho,obj.wt,obj.maxDeltaX);
            thetaMin = sqrt(2)*erfinv(MCS.F)*theta0;
            t = (2*obj.particle.momentum*sin(thetaMin/2))^2;
        end
        
        function [dE,dAngles,newSecondaryInfo] = interact(obj)
            % sampling procedure
            tm = obj.tmin;
            accepted = false;
            while ~accepted
                z = rand(2,1);
                t = tm/(1-z(1));
                % accept t with probability g(t)
                a = 0.856e3*obj.Rsquared;
                gNormalised = exp(-a*(t-tm));
                if z(2)<gNormalised
                    accepted = true;
                    t = sym(t);
                end
            end
            mn = sym(obj.MN);
            mp = sym(obj.particle.mass);
            pp1 = sym(obj.particle.momentum);
            Ep1 = obj.particle.energy;
            En2 = 0.5*t/mn + mn;
            pn2 = sqrt(En2^2-mn^2);
            Ep2 = Ep1-t/(2*mn);
            pp2 = sqrt(Ep2^2-mp^2);
            cosTheta = (pp1^2+pp2^2-pn2^2)/(2*pp1*pp2);
            sinTheta = sin(acos(cosTheta));
            phi = 2*pi*rand;
            dAngles = [sin(phi);cos(phi);double(sinTheta);double(cosTheta)];
            dE = double(Ep2-Ep1);
%             obj.particle.energy = double(Ep2);
            newSecondaryInfo = double.empty(6,0);
        end

% function [dE,dAngles,newSecondaryInfo] = interact(obj)
%             tm = obj.tmin;
%             accept = false;
%             gmin = obj.gMott(tm);
%             gmax = obj.gMott(Inf);
%             while ~accept
%                 z = rand(2,1);
%                 tdouble = tm./(1-z(1)); % draw from pdf 1/t^2 from tmin to inf
%                 g = (obj.gMott(tdouble)-gmin)/(gmax-gmin);
%                 if z(2)<g
%                     accept = true;
%                     t = sym(tdouble);
%                 end
%             end
%             
%         
%         end
    end

    methods(Access=private)
        function gt = gMott(obj,t)
            gamma = obj.particle.gamma;
            beta = obj.particle.beta;
            gt = 1-beta^2*t/(4*(Consts.megev*beta*gamma));
        end
    end
end

