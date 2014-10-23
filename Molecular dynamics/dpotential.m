%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function of derivative of potential
% given by -12[1/r^13-1/r^7]

function dp = dpotential(r)

dp = -12*(1./r.^13-1./r.^7);

end