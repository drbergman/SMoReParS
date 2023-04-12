function out = computeTimeSeries(p,tt,~,~)

opts = odeset("NonNegative",1);
[~,out] = ode23(@(t,x) odefn(x,p),tt,100,opts);

out(isnan(out)) = -1e6; % something really different from whatever data I'm comparing it to so that it gets heavily penalized for getting NaN

if length(out)<length(tt) % for some parameter values, the solution blows up in finite time and this penalized those
    out((length(out)+1):length(tt)) = -1e6;
end
