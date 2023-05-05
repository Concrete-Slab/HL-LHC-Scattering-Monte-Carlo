classdef NuclearProcess < HardProcess
    %NuclearProcess Summary of this class goes here
    %   Detailed explanation goes here

    properties
        A(1,1) double       % average relative atomic mass
        neff(1,1) double    % effective nucleon number for incoherent scattering

    end

    properties(Dependent)
        sigmaTot
    end

    properties(Dependent,SetAccess=private)
        sigmapnElastic
        sigmapnSD
        sigmapnTotal
        sigmapNTotal
        sigmapNElastic
        sigmapNInelastic
        bpp % slope parameter for elastic proton-proton differential sigma
    end

    methods
        function obj = NuclearProcess(material,particle)
            arguments
                material(1,1) Material = Material
                particle(1,:) ParticleHandle {mustBeScalarOrEmpty} = ParticleHandle.empty
            end
            %NuclearProcess Construct an instance of this class
            %   Detailed explanation goes here
            obj@HardProcess(material,particle);
            obj.A = sum(material.wt.*material.A);
            obj.neff = 1.1618*obj.A^(1/3);
        end

        function st = get.sigmaTot(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            if obj.particle.mass == Consts.mpgev && obj.particle.charge == 1
                st = obj.sigmapnSD + obj.sigmapnElastic + obj.sigmapNElastic;
                % convert mbarn to m^2
                st = obj.sigmapNTotal;
            else
                st = 0;
            end
        end

        function [dE,dAngle,newSecondaryInfo] = interact(obj)

            newSecondaryInfo = double.empty(6,0);
            dE = 0;
            dAngle = 0;
            % sample the nuclear process that happens from their relative
            % cross-sections
            % no differential data for inelastic - it just adds to total
            % cross-section
            cmf = [0,obj.sigmapnElastic,obj.sigmapnSD,obj.sigmapNElastic];
            for i = 2:length(cmf)
                cmf(i) = cmf(i)+cmf(i-1);
            end
            cmf = cmf./cmf(end);
            if any(~isfinite(cmf)) || any(~isreal(cmf))
                dE = 0;
                dAngle = [0;1;0;1];
                return
            end
            zProcess = rand;
            I = find(zProcess<cmf,1,"first");
            p1 = obj.particle.momentum;
            E1 = obj.particle.energy;
            fmInitLab = [E1;0;0;p1]; % z directional motion
            switch I
                case 2
                    % proton-nucleon elastic process
                    t = sampleppElastic(obj);
                    [betaCOM,~] = getCOM(p1,E1,0,Consts.mpgev);
                case 3
                    % proton-nucleon single diffractive scattering
                    [dE,dAngle] = samplepnSD(obj);
                    return
                case 4
                    t = samplepNElastic(obj);
                    [betaCOM,~] = getCOM(p1,E1,0,obj.A*Consts.mpgev); % approximate nucl. mass as being A protons
                otherwise
                    dE = 0;
                    dAngle = [0;1;0;1];
                    return
            end
            % perform frame changes for elastic process
            fmInitCOM = lorentz(fmInitLab,betaCOM);
            pCOMMag = norm(fmInitCOM(2:4));
            ECOM = fmInitCOM(1);
            thetaCOM = sym(sqrt(t/pCOMMag));
            phi = 2*pi*rand;
            sinPhi = sym(sin(phi));
            cosPhi = sym(cos(phi));
            pFinalCOM = pCOMMag*double([cosPhi*sin(thetaCOM);sinPhi*sin(thetaCOM);cos(thetaCOM)]);
            fmFinalCOM = [ECOM;pFinalCOM];
            fmFinalLab = inverselorentz(fmFinalCOM,betaCOM);
            dE = fmFinalLab(1)-fmInitLab(1);
            symbolicMomentum = sym(fmFinalLab(2:4)); % use variable precision arithmetic
            theta = acos(symbolicMomentum(3)/norm(symbolicMomentum));
            dAngle = double([sinPhi;cosPhi;sin(theta);cos(theta)]); % convert to double
        end
    end

    methods % getters, private
        function sppel = get.sigmapnElastic(obj)
            p = obj.particle.momentum;
            sppel = obj.neff*7.0e-31*(p/450)^(4.79e-2); % barn
        end

        function sppsd = get.sigmapnSD(obj)
            sppsd = obj.neff*0.68e-31*log(0.15*2*Consts.mpgev*obj.particle.momentum); % barn

        end

        function spptot = get.sigmapnTotal(obj)
            p = obj.particle.momentum;
            spptot = 40e-31*(p/450)^(5.79e-2);
        end

        function spNtot = get.sigmapNTotal(obj)
            sigma0T = 50.59e-31;
            alpha = 0.77;
            spNtot = sigma0T*obj.A^alpha;
        end

        function spNel = get.sigmapNElastic(obj)
            % elastic pN sigma is total minus every other sigma
            spNel = obj.sigmapNTotal-obj.sigmapNInelastic-obj.sigmapnTotal;

        end

        function spNinel = get.sigmapNInelastic(obj)
            sigma0I = 41.2e-31;
            alpha = 0.711;
            spNinel = sigma0I*obj.A^alpha;
            % correlation assumes energy independence at high energy, data from:
            % https://deepblue.lib.umich.edu/bitstream/handle/2027.42/23455/0000406.pdf%3Bjsessionid%3D5141B4FA6A8B33F128B3BBDDF396D83D?sequence%3D1
        end

        function b = get.bpp(obj)
            s = 2*Consts.mpgev*obj.particle.momentum;
            b = 8.5+1.086*log(sqrt(s));
        end
    end

    methods % sampling distributions
        function [t] = sampleppElastic(obj)
            % centre of momentum energy - Minkowski variable
            %             p = obj.particle.mass;
            %             E = obj.particle.energy;
            %             s = 2*Consts.mpgev*p;
            %             bpp = 8.5+1.086*log(sqrt(s));
            b = obj.bpp;
            tmax = 1.2e-2;
            z = rand;
            t = -1/b * log(1-z*(1-exp(-b*tmax))); % sample t from exponential dist
            %             [betaCOM,~] = getCOM(p,E,0,Consts.mpgev);
            %             fmInitialLab = [E;0;0;p];
            %             fmInitialCOM = lorentz(fmInitialLab,betaCOM);
            %             momCOMmag = norm(fmInitialCOM(2:4));
            %             ECOM = fmInitialCOM(1);
            %             thetaCOM = sqrt(t/momCOMmag); % t = (p*theta)^2
            %             phi = 2*pi*rand;
            %             pFinalCOM = momCOMmag*[cos(phi)*sin(thetaCOM);sin(phi)*sin(thetaCOM);cos(thetaCOM)];
            %             fmFinalCOM = [ECOM;pFinalCOM];
            %             fmFinalLab = inverselorentz(fmFinalCOM,betaCOM);
            %             momFinalLab = sym(fmFinalLab(2:4));
            %             EFinalLab = sym(fmFinalLab(1));
            %             theta = acos(momFinalLab(3)./norm(momFinalLab));
            %             dAngles = [sin(phi),cos(phi),double(sin(theta)),double(cos(theta))];
            %             dE = double(EFinalLab-sym(E));
        end

        function [dE,dAngles] = samplepnSD(obj)
            % first sample M:
            EiLab = obj.particle.energy;
            s = sym(2*Consts.mpgev*EiLab+2*Consts.mp^2);
            zM = rand;
            % approximation to marginal for Msquare is p(M^2)=1/M^2.
            % the cited paper in the report seems to have assumed 0<t<inf
            % when calculating the marginal for M, so I have done the same
            Msquare = (0.15*s)^zM;
            % now sample t given Msquare. p(t|M) = b*exp(b*(M^2m_p^2/s-t))
            % again use inverse transform sampling
            b = 7/12 * obj.bpp;
            zt = rand;
            tmax = 0.1;
            % find tmin from CM considerations
            M = sqrt(Msquare);
            % energy change from releasing mass M in CM
            delE = sym((M^2+2*M*Consts.mpgev)/(2*sqrt(s)));
            % initial momentum in CM
            piCM = vpa(sqrt(s/4 - Consts.mpgev^2));
            % final momentum in CM
            pfCM = sqrt((sqrt(s)/2 - delE)^2-Consts.mpgev^2);
            % minimum t is determined by kinematics of CM
            tmin = delE^2-(piCM^2+pfCM^2+2*piCM*pfCM);
            % sample t from its distribution
            t = -1/b * log(exp(-b*tmin)-zt*(exp(-b*tmin)-exp(-b*tmax)));
            % use kinematics to get cosTheta in CM frame from t
            cosThetaCM = 1-(t-tmin)/(2*piCM*pfCM);
            sinThetaCM = sin(acos(cosThetaCM));
            % random azimuth, doesnt change with z lorentz boost
            phi = rand*2*pi;
            % perform inverse lorentz boost back to lab frame
            fMomentumVecCM = pfCM.*[cos(phi)*sinThetaCM;sin(phi*sinThetaCM);cosThetaCM];
            fmfCM = [sqrt(pfCM^2+Consts.mpgev^2);fMomentumVecCM];
            betaCM = getCOM(obj.particle.momentum,EiLab,0,Consts.mpgev);
            fmfLab = sym(inverselorentz(double(fmfCM),betaCM));
            % obtain final results
            cosThetaLab = fmfLab(4)./norm(fmfLab(2:4));
            sinThetaLab = sin(acos(cosThetaLab));
            dE = double(EiLab-fmfLab(1));
            dAngles = double([sin(phi);cos(phi);sinThetaLab;cosThetaLab]);
        end

        function t = samplepNElastic(obj)
            bpN = 14.1*obj.A^0.65; % GeV^-2
            tmax = 1.2e-2;
            z = rand;
            t = -1/bpN * log(1-z*(1-exp(-bpN*tmax)));
        end
    end
end

