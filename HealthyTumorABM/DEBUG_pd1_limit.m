function DEBUG_pd1_limit(M)

assert(all(M.checkpoint.PD1 <= M.checkpoint.pd1_on_immune * (1+eps)));
assert(all(M.checkpoint.PD1 + M.checkpoint.PD1aPD1 <= M.checkpoint.pd1_on_immune * (1.0001)));
