function [kbounds] = pk_bounds(u,M,K,Qip,neig)
% [kbounds] = pk_bounds(M,K,Qip,ub,qoo,neig)
%
% Compute upper/lower frequency bounds for the pk-bisection solver.
%
% Mr:        = (u/b)^2 * M : reduced mass matrix
% K:         stiffness matrix
% Qip:       aerodynamic load interpolation table (struct)
% neig:      number of aeroelastic eigenvalues to use
% kbounds:   frequency limit values for bisection (vector)
%
% (c) 2016 Dan Borglund <dodde@kth.se>, Ulf Ringertz <rzu@kth.se>
%          and David Eller <dlr@kth.se>

  % reduced frequency to start from (small, but > 0)
  kinit = 0.01;

  % check if options make sense
  sdim = size(M, 1);
  if neig > sdim
    error('pk_bounds: Requested neig > system dimension.');
  end

  kbounds = zeros(1, neig+1);   % for storage of bounds
  kbounds(1) = kinit;           % first lower bound if kinit is small
  kstep = kinit;                % initial step in k
  k = kinit-kstep;              % initiate k
  nbound = 0;                   % number of bounded eigenvalues

  % reduced mass matrix
  b = Qip.bref;
  Mr = (u/b)^2 * M;

  % dynamic pressure
  qoo = 0.5*Qip.rho*u^2;

  while nbound < neig

    k = k + kstep;               % update k

    % retrieve Q(k)
    Qk = ipolQk(Qip, k);

    % solve p-k eigenvalue problem
    pk2 = eig(qoo*Qk - K, Mr);        % solve ((u2/b2)*M*p2+K-qA(k))v=0
    pk = sqrt(pk2);                % solve for p
    pk = sign(imag(pk)).*pk;       % we want the square roots p with k>0
    impk = sort(imag(pk));         % sort p in frequency

    nbelow = sum(impk < k);        % number of eigenvalues with im(p)<k
    if nbelow == nbound            % no crossing of im(p)=k has occured
      kstep = 2*kstep;             % increase step in k
    elseif nbelow == (nbound + 1)  % one crossing has occured
      nbound = nbound+1;           % update number of found bounds
      kbounds(nbound+1) = k;       % store bound and keep kstep
    else                           % more than one crossing has occured
      k = k-kstep;                 % start over
      kstep = kstep/2;             % with smaller step in k
    end

%     % check if k is sane
%     if k > 100.0
%       error('pk_bounds failed: unreasonable reduced frequency.');
%     end
  end
