function [delta,ec50] = computeHillParameters(d,delta_d,b)

if d(1)==0
    d(1) = .1*d(2);
end
C = [.75;7.55];
D = d + delta_d;

b_log_ec50 = log(delta_d) - log(d./C(1)^b-D/C(2)^b);
log_ec50 = b_log_ec50/b;
ec50 = exp(log_ec50);
delta = d.*((ec50/C(1)).^b+1);

if ~isreal(delta) || ~isreal(ec50)

    error("Complex parameters.")

end