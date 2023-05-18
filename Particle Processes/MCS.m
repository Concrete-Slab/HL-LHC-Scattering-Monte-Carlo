classdef MCS < SoftProcess
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(SetAccess=private)
        X0(1,:) double          % radiation length (g cm^-2)
        A(1,:) double           % atomic mass (AMU)
        Z(1,:) double           % atomic number
        wt(1,:) double          % mass percentage
        rhoi(1,:) double        % density of each component (g cm^-3)
    end

    properties(Access=private)
        fit(1,1) string {mustBeMember(fit,["LynchDahl","Highland"])} = "LynchDahl"
    end

    properties(Dependent,SetAccess=protected)
        dEdx
    end

    properties(Constant)
        s2(1,1) double = 13.6
        epsilon(1,1) double = 0.088;
        F(1,1) double = 0.98
    end
    
    methods
        function obj = MCS(material)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            arguments(Input)
                material(1,1) Material;
            end
            obj@SoftProcess
            obj.X0 = material.X0*0.1;
            obj.A = material.A;
            obj.Z = material.Z;
            obj.rhoi = 0.001*material.rho*material.wt;
            obj.fit = Consts.data.mcsFit;
        end
        
        function [delE,delAngles,delPos] = update(obj,dr)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            arguments(Input)
                obj(1,1) MCS
                dr(1,1) double {mustBeNonnegative}  % path length, m
            end
            arguments(Output)
                delE(1,1) double                    % energy change in GeV/c
                delAngles(2,1) double               % phi,theta in local frame
                delPos(2,1) double                  % change in transverse position, m
            end
            % Elastic collision
            delE = 0;
            % get particle motion values
            p = obj.particle.momentum*1000; % convert to MeV/c in natural units from GeV/c in natural units
            zp = obj.particle.charge;
            beta = obj.particle.beta;
            % determine std. deviation using the chosen fit
            if obj.fit == "Highland"
                % Convert path length into g cm^2
                X = sum(obj.rhoi)*(100*dr);
                % apply highland correlation
                sigma = 13.6./(p*beta)*sqrt(X./obj.X0)*(1+0.088*log10(X./obj.X0));

            else

                % path length in g/cm^2 for each component
                Xi = obj.rhoi.*(100*dr);
                
                chiAi = 2.007e-5.*obj.Z.^(2/3).*(1+3.34.*(obj.Z.*zp.*Consts.alpha./beta).^2)./(p^2);
                chiCi = 0.157*(obj.Z.*(obj.Z+1).*Xi./obj.A).*((zp/(beta*p))^2);
                
                chiA = exp(sum((Xi.*obj.Z.*(obj.Z+1)./obj.A).*log(chiAi))./sum(Xi.*obj.Z.*(obj.Z+1)./obj.A));
                chiC = sum(chiCi);

                omega = chiC/chiA;
                v = 0.5*omega/(1-obj.F);
                sigma = sqrt(chiC/(1+obj.F^2)*((1+v)/v * log(1+v) -1));
            end
            % Generate random normal variables
            zn = randn([1,4]);
            % apply correlations
            thetax = zn(1)*sigma;
            thetay = zn(2)*sigma;
            dx = zn(3)*dr*sigma./sqrt(12) + dr*thetax./2;
            dy = zn(4)*dr*sigma./sqrt(12) + dr*thetay./2;
            % set outputs
            delPos = [dx;dy];



            % convert from plane angles to polar angles
            if thetay==0
                phi = 0;
            else
                phi = sign(tan(thetay))*acos(tan(thetax)/sqrt(tan(thetax)^2+tan(thetay)^2));
            end
%             theta = acos(1./sqrt(1+tan(thetax).^2+tan(thetay).^2));
            % appropriate if -pi<theta<pi
            theta = atan(sqrt(tan(thetax).^2+tan(thetay).^2));
            
            

            delAngles = [phi;theta];
            
%             delAngles = rand(2,1);

%             delAngles = [thetax;thetay];
        end

        function e = get.dEdx(obj)
            e=0;
        end
    end

    methods(Static)
        function zNormal = sampleCentralF(theta0,n)
            %TODO implement this function
            z = rand(1,n);
            thetaStandard = sqrt(2)*erfinv(MCS.F*(2*z-1));
            % scale by the standard deviation
            zNormal = theta0*thetaStandard;
        end

        function thetaMax = getThetaMax(theta0)
            thetaMax = sqrt(2)*erfinv(MCS.F)*theta0;
        end

        function theta0 = getTheta0(particleMomentum,particleBeta,particleCharge,Zi,Ai,X0,rho,wt,deltaX)
            arguments
                particleMomentum(1,1) double    % GeV/c
                particleBeta(1,1) double        % 1
                particleCharge(1,1) double      % e
                Zi(1,:) double                  % e
                Ai(1,:) double                  % AMU
                X0(1,1) double                  % kg m^-2
                rho(1,1) double                 % kg m^-3
                wt(1,:) double                  % 1
                deltaX(1,1) double              % m
            end
            
            rhoi = 0.001*rho*wt; % convert kg/m^3 to g/cm^3
            deltaX = 100*deltaX; % convert m to cm
            X0 = 0.1*X0; % convert kg/m^2 to g/cm^2
            fit = Consts.data.mcsFit;
            p = particleMomentum*1000; % convert GeV/c to MeV/c;
            if fit == "Highland"
                % Convert path length into g cm^2
                X = sum(rhoi)*(deltaX);
                % apply highland correlation
                theta0 = 13.6./(p*particleBeta)*sqrt(X./X0)*(1+0.088*log10(X./X0));
                if ~isreal(theta0)
                    theta0
                end
            else

                % path length in g/cm^2 for each component
                Xi = rhoi.*(100*deltaX);
                
                chiAi = 2.007e-5.*Zi.^(2/3).*(1+3.34.*(Zi.*particleCharge.*Consts.alpha./particleBeta).^2)./(p^2);
                chiCi = 0.157*(Zi.*(Zi+1).*Xi./Ai).*((particleCharge/(particleBeta*p))^2);
                
                chiA = exp(sum((Xi.*Zi.*(Zi+1)./Ai).*log(chiAi))./sum(Xi.*Zi.*(Zi+1)./Ai));
                chiC = sum(chiCi);

                omega = chiC/chiA;
                v = 0.5*omega/(1-MCS.F);
                theta0 = sqrt(chiC/(1+MCS.F^2)*((1+v)/v * log(1+v) -1));
            end
        end
    end
end

