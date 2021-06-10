function [urev,zrev] = reversal(K,Qip,f,CRv,CRd)
% [urev,zrev] = reversal(dsy, geo)
%
% Compute the airpeed at aileron reversal and the deformation
% for unit flap deflection at reversal speed.
%
% geo:    wing geometry description
% dsy:    aeroleastic system
% urev:   critical airspeed at aileron reversal
% zrev:   wing deformation vector for unit delta at qrev
%
% (c) 2016 Dan Borglund <dodde@kth.se> and David Eller <dlr@kth.se>

% Qip   : struct containing Q(k) precomputed at a set of reduced frequencies
% type 'help ipolQk' at the prompt for an explanation of the format

% f     : aileron load vector per unit deflection
% CRv   : coefficient of rolling moment due to displacement
% CRd   : coefficient of rolling moment due to aileron deflection

Q_delta = -f*CRv/CRd;
Q_0 = Qip.Qtab;
Q_0 = Q_0(:,:,1);
Q = Q_0 + Q_delta;
[V,D] = eig(Q,K);
% ev = abs(diag(D));
ev = diag(D);
ii = find(imag(ev) ~= 0);
ev(ii) = 0;
qrev = 1/max(ev);
index = find(max(ev) == ev);
rho = Qip.rho;

urev = sqrt(2*qrev/rho);

zrev = V(:,index);

end

% should I see twist when I plot the modeshape?

% the reversal speed is lower wrt the divergence one, I think it makes
% sense, otherwise we would only experience divergence then the wing will
% break (?)
