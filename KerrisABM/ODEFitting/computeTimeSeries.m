function out = computeTimeSeries(p,tt,~,~)

if p(3)>p(1) % beta should not be bigger than alpha (or else the solution will immediately start to decrease for any x>1 since theta<1)
    out = 100*exp((p(1)-p(3))*tt'); % just approximate this as exponential decay to avoid issues of the solution going negative and giving complex answers that even NonNegativity constraints cannot fix
    return;
end
p(2) = 1 - 1/p(2); % p(2) is nu coming in and theta going forward
opts = odeset("NonNegative",1);
[~,out] = ode45(@(t,x) odefn(x,p),tt,100,opts);

if any(isnan(out))
    out(isnan(out)) = -1e6; % something really different from whatever data I'm comparing it to so that it gets heavily penalized for getting NaN
end

if length(out)<length(tt) % for some parameter values, the solution blows up in finite time and this penalized those
    out((length(out)+1):length(tt)) = -1e6;
end

out = min(out,1e10); % for plenty of chosen parameters (alpha >> beta for example) the solutions can grow large, but finite, s.t. the error they produce is Inf
