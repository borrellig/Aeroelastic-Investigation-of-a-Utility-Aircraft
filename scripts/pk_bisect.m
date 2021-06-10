function [pu,vu] = pk_bisect(ieig,u,M,K,Qip,kbounds)
% [pu,vu] = pk_bisect(u,ieig,kbounds,dsy,geo,opt)
%
% Solve the nonlinear pk-eigenvalue problem using bisection.
% ieig is the index of the eigenmode considered, in order of
% increasing frequency (ieig = 1 corresponds to the eigenmode
% with lowest reduced frequency k).
% This functions requires the option opt.ktol, which is the
% required accuracy for the reduced frequency.
%
% u:         airspeed [m/s]
% ieig:      eigenvalue index (integer)
% kbounds:   frequency limits as obtained from pk_bounds (vector)
% dsy:       aeroelastic system (struct)
% geo:       wing description (struct)
% opt:       solver options (struct)
% vu:        aeroelastic eigenvector at u (vector)
% pu:        flutter eigenvalue at u (scalar)
%
% (c) 2016 Dan Borglund <dodde@kth.se>, Ulf Ringertz <rzu@kth.se>
%          and David Eller <dlr@kth.se>
  
  % tolreance for iterative solution of pk-eigenvalue problem
  ktol = 0.001;

  % reduced mass matrix
  b = Qip.bref; % aeroelastic reference length (half chord)
  Mr = (u/b)^2 * M;

  % dynamic pressure
  qoo = 0.5*Qip.rho*u^2;

  % upper/lower bounds
  klow = kbounds(ieig);
  kupp = kbounds(ieig+1);

  converged = 0;
  while not(converged)

    % new estimate for k
    k = (klow+kupp)/2;

    % Fill in code to solve the flutter eigenvalue problem at the
    % current approximation for the reduced frequency k.
    % Your code must compute the scalar variable p corresponding
    % to the eigenvalue with index ieig (argument) and the corresponding
    % complex eigenvector v.
    
    Q = ipolQk(Qip,k);
    [V,D] = eig(K-qoo*Q,-Mr);
    p2 = diag(D);
    p = sqrt(p2);
    p = p.*sign(imag(p));
    [B,I] = sort(imag(p));
    index = I(ieig);
    p = p(index);
    v = V(:,index);
    
    imp = imag(p);                 % imaginary part of current p
    if imp > k                     % then k new lower bound
      klow = k;                    % update lower bound
    elseif imp <= k                % then k new upper bound
      kupp = k;                    % update upper bound
    end
    converged = (kupp-klow) < ktol;  % eigenvalue bounded to within ktol?
  end
  
  
  % fill return values
  pu = p;
  vu = v;
