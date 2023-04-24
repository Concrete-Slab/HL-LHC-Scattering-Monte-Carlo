function dEdx = betheMatt(M, A, Z, I, w, rho, beta, z, deltaOn)
% BETHE Bethe equation - Matt Lenahan 2023
%   A           1*n vector of atomic mass [g mol-1]
%   Z           1*n vector of atomic no.
%   I           1*n vector of mean ionisation energy [MeV]
%   w           1*n vector of mass fractions
%   rho         1*n vector of density [g cm-3]
%   beta        v/c
%
%   dEdx        mass stopping power [MeV g cm^-2]

    % CONSTANTS
    K = 0.307075; % MeV mol-1 cm2
    m_e = .51099895000; % electron mass MeV/c^2
    %z = 1; % incident proton charge
    %M = 938.27208816; % incident proton mass MeV/c^2

    gamma = 1/sqrt(1-beta^2);
    W_max = wmax(beta, gamma, M); % MeV
    I = bragg(A, Z, I, w); % MeV
    if deltaOn == 1
        delta = deltacorr2(beta*gamma, w, Z, A, I, rho);
    else
        delta = 0;
    end

    ZoverA = sum(w.*Z./A); % see D.E. Groom et al. (2001) app. A
    
    t1 = K*z^2*ZoverA*beta^-2;
    t2 = .5*log(2*m_e*beta^2*gamma^2*W_max*I^-2);
    
    dEdx = t1*(t2-beta^2-delta/2);

end

function W_max = wmax(beta, gamma, M)
% WMAX

    % CONSTANTS
    m_e = .51099895000; % electron mass MeV/c^2
    %M = 938.27208816; % incident proton mass MeV/c^2
    
    t1 = 2*m_e*beta^2*gamma^2;
    t2 = 1+(2*gamma*(m_e/M))+(m_e/M)^2;
    W_max = t1/t2;
end

function I = bragg(A, Z, I, w)
%BRAGG calculate I for compounds

    n = sum(w.*(Z./A).*log(I));
    d = sum(w.*(Z./A));
    I = exp(n/d);

end

function delta = deltacorr2(betagamma, w, Z, A, I, rho)
%DELTACORR2 delta correlation based on D.E. Groom et al. (2001)

    I = I*1e6; % MeV to eV
    ZoverA = sum(w.*Z./A); % see appendix A
    hwp = 28.816*(rho*ZoverA)^.5; % eV
    x = log10(betagamma);
    C = 2*log(I/hwp) + 1;
    if I < 100
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
    
    t1 = C - 2*log(10)*x0;
    t2 = (x1-x0)^3;
    a = t1/t2;

    if x >= x1
        delta = 2*log(10)*x - C;
    elseif x >= x0
        delta = 2*log(10)*x - C + a*(x1-x)^3;
    else
        delta = 0;
    end

end

