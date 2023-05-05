function f = normalDistribution(mu,sigma)
    f = @(x) (1/sigma^2*2*pi).*exp(-0.5*((x-mu)./sigma).^2);
end

