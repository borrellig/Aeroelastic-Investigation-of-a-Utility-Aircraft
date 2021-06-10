function Qk = ipolQk(Qip, k)
% Qk = ipolQk(Qip, k)
%
% Perform piecewise linear interpolation of the aerodynamic loads
% stored in Qip, which is a struct with (at least these) two elements:
%
% Qip.ktab : is a vector of the nk reduced frequencies at which Q(k) has been
%            precomputed and stored in
% Qip.Qtab : a 3-dimensional array with size [n,n,nk] which contains nk
%            n-by-n load matrices (last index is frequency index)
% Qk : the interpolated load matrix Q(k)
%
% (c) 2004-2016 Dan Borglund <dodde@kth.se> and David Eller <dlr@kth.se>

%   if k < 0
%     Qk = conj(ipolQk(Qip, -k));
%     return;
%   end

  kt = Qip.ktab;
  Qt = Qip.Qtab;
  nk = length(kt);
  ndof = size(Qt,1);

  % cap Q(k) above and below min and max(ktab), do not extrapolate
  if k <= kt(1)
    Qk = reshape(Qt(:,:,1), [ndof, ndof]);
  elseif k >= kt(nk)
    Qk = reshape(Qt(:,:,nk), [ndof, ndof]);
  else

    % search for k in ktab
    ilo = 1;
    ihi = nk;
    while ihi-ilo > 1
      imid = fix((ilo+ihi)/2);
      kmid = kt(imid);
      if kmid <= k
        ilo = imid;
      else
        ihi = imid;
      end
    end

    % linear interpolation
    t = (k - kt(ilo)) / (kt(ihi) - kt(ilo));
    Qlo = reshape(Qt(:,:,ilo), [ndof, ndof]);
    Qhi = reshape(Qt(:,:,ihi), [ndof, ndof]);
    Qk = (1-t)*Qlo + t*Qhi;

  end

end
