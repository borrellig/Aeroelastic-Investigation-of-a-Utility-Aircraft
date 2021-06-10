function [udiv,zdiv] = divergence(K,Qip)
% [udiv,zdiv] = divergence(K,Qip)
%
% Compute the divergence speed and the corresponding
% deformation shape.
%
% K:     stiffness matrix
% Qip:   aerodynamic load interpolation struct
% udiv:  divergence speed
% zdiv:  divergence deformation shape
%
% (c) 2004-2016 Dan Borglund <dodde@kth.se> and David Eller <dlr@kth.se>

  % solve numerically well-posed divergence eigenvalue problem
  % for the critical dynamic pressure qdiv

  % for the steady state analysis we consider the case when
  % k = 0 (reduced frequency) so to have no flapping. First option
  % (Qip.ktab(1))
  A = Qip.Qtab;
  A = A(:,:,1);
  [V,D] = eig(A,K);
  ev = diag(D);
  ii = find(imag(ev) ~= 0); % we only need the max between real eigenvalues
  ev(ii) = 0;
  
  
  % ev = diag(D)
  qd = 1/max(ev);
  index = find(max(ev) == ev);
  rho = Qip.rho; % air density
  
  % determine critical speed from dynamic pressure
  udiv = sqrt(2*qd/rho);

  % corresponding eigenvector
  % zdiv = V(:,index);
  zdiv = V(:,index);

end

% !!! remember !!!
% in order to have divergence the matrix K has to be non singxular, thus we
% have to add the boundary conditions for a clamped side, for instance,
% otherwise udiv = 0;
