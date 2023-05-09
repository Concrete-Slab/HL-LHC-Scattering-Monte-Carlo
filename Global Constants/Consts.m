classdef Consts
    %CONSTS Static constants for use by all scripts
    %   Global constants used by the simulator are stored as properties.
    %   This saves use of the workspace, which is volatile, the need to
    %   save constants in multiple different files, or additional
    %   unneccessary function arguments
    
    properties(Constant)
        % minkowski norm
        mNorm(4,4) double = [1 0 0 0;0 -1 0 0;0 0 -1 0;0 0 0 -1]
        % proton rest mass (kg)
        mp(1,1) double = 1.67262192369e-27
        % proton rest mass (GeV/c^2)
        mpgev(1,1) double = Consts.c^2*Consts.mp/(1e9*Consts.e);
        % electron rest mass (kg)
        me(1,1) double = 9.109383701528e-31
        % electron rest mass (GeV/c^2)
        megev(1,1) double = Consts.c^2*Consts.me/(1e9*Consts.e);
        % speed of light in a vacuum (m/s)
        c(1,1) double = 299792458
        % elementary charge (C)
        e(1,1) double = 1.602176634e-19
        % incident proton four-momentum (GeV/c), natural units
        initialFourMomentum(4,1) double = [7000 0 0 6.999999937117533e+03]
        % Avogadros constant, /kmol
        Na = 6.02214076e26;
        % vacuum permittivity
        e0 = 8.854187812813e-12;
        % classical electron radius, m
        re = Consts.e^2./(4*pi*Consts.e0*Consts.me*Consts.c^2)
        % fine structure constant
        alpha = 1/137;
        % Molar Mass constant kg kmol^-1
        MolarMassConstant =  0.9999999996530
        % reduced planck constant * speed of light (J m)
        hbarc = 3.16152677e-26

        % Bethe-Bloch constant K (GeV kmol^-1 m^2)
        Kbb = 4*pi*Consts.Na*Consts.re^2*Consts.megev
        
        % LHC Constants
        % LHC magnet material
        MagnetMaterial = "NbTi"
        % LHC beam tube diameter is about 0.025m
        TubeRadius = 0.025;
        % maximum beam sigma and mu
        MaxSigma = 240e-6;
        MaxMu = 0.002;
        TungstenDepth = 0.016;
        % data handle, values are global and variable
        data(1,1) RuntimeData = RuntimeData;
    end

    properties(Constant,Access=private)
        % minimum kinetic energies for protons and electrons
        MKEElectron = 1e-3; % (25e-6 GeV = 1MeV) limit of D.Tan relation
        MKEProton = 10e-6; % (10e-6 GeV = 10keV) from geant4 suggestion
    end

    methods(Static)
        % minimum kinetic energy to propagate a particle, GeV
        function mke = MinimumKineticEnergy(mass,charge)
            arguments
                mass(1,:) double
                charge(1,:) double 
            end
            protonIdx = and(mass==Consts.mpgev,charge==1);
            electronIdx = and(mass == Consts.megev,charge ==-1);
            mke = zeros(1,length(mass));
            mke(protonIdx) = Consts.MKEProton;
            mke(electronIdx) = Consts.MKEElectron;
            if any(mke==0)
                warning("Unrecognised particle")
            end
            
        end
    end
end

