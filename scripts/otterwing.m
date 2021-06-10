function [M,K,Z,Qip,f,CRv,CRd,EI,GK] = otterwing(B, l, b, t, ba, mhinge, rhop, E, G, nelem, dpm, Ac,hspar)
% [M,K,Z,Qip,f,CRv,CRd] = labwing(B, l, b, t, ba, mhinge, rhop, E, G, nelem, dpm)
%
% Assemble structural system matrices and aerodynamic load interpolation.
%
% Input parameter:
% l     : semispan [m]
% b     : semichord [m]
% ba    : aileron semichord [m]
% mhinge: mass of hinge axis [kg]
% t     : plate thickness [m]
% rhop  : plate density [kg/m^3]
% E     : Flexural modulus of elasticity [N/m^2]
% G     : Shear modulus [N/m^2]
% nelem : number of finite elements
%
% dpm   : Matrix used to define a set of concentrated masses to attach to the
%         structural model. dpm has 3 columns and one row for each concentrated mass
%         such that each row contains mass, x-position and y-position:
%         [ m1 xpos1 ypos1
%           m2 xpos2 ypos2 ...
% B     : Matrix of homgeneous linear constraints, one constraint per row.
%
% Outputs:
% M     : mass matrix
% K     : stiffness matrix
% Z     : null-space basis for the constraints B. M, K, Qip are already projected.
% Qip   : struct containing Q(k) precomputed at a set of reduced frequencies
%         type 'help ipolQk' at the prompt for an explanation of the format
% f     : aileron load vector per unit deflection
% CRv   : coefficient of rolling moment due to displacement
% CRd   : coefficient of rolling moment due to aileron deflection
%
% (c) 2004-2016 Dan Borglund <dodde@kth.se> and David Eller <dlr@kth.se>

  % at least, warn about inappropriate discretization
  if mod(nelem,3) ~= 0
    fprintf(1,'Warning: Number of finite elements should be divisible by 3.\n');
  end

  % define properties for 4 spanwise support points
  % nsup = 4;
  nsup = nelem+1;
  % y = linspace(0,l,nsup); % !!! if test with clamped
  y = l/2*linspace(-1, 1, nsup); % so to allow for negative y-coord when 
  % positioning concentrated masses in dpm

  % planform definition
  chord = 2*b;
  chail = 2*ba;
  xle = -b;
  xte = +b;
  ca = chail*ones(nsup,1); % accounts for ailerons
  ca(0.25*nsup+1:0.75*nsup) = 0;
  % only on outer half of the wing
  % ca = chail*ones(nsup,1);
  % xea = -0.5*chail;
  xea = -0.15*chord; % el. axis at 35% of chord (coordinate wrt half chord)

%   % mass and inertia properties
%   mp = rhop*chord*t;     % wing plate mass/unit span [kg/m]
%   xp = 0;                % plate chordwise center of mass [m]
%   Jp = (mp/12)*chord^2;  % plate rotary inertia wrt plate cm [kgm] %
%   small thickness approximation (compare with data sheets)

  % wingbox properties, NACA 63A-516
  [l_wbox,xle_wbx,xte_wbx,le_wbx,te_wbx,x_ud,l_up,l_down] = wbox_prop(b);
        
  % mass and inertia properties of wing
  A_wbox = t*l_wbox;
  m_wbox = A_wbox*rhop;
  
  % tapered spar
  l_e = l/(nelem);
  tspar = 0.003; % change it to 0.005 for constant Ac
  mspar = (sum(Ac*l_e)/l+0.14*chord*tspar)*rhop;
  msparfun = ((Ac*l_e)/l+0.14*chord*tspar)*rhop;
  mp = mspar+m_wbox;
  % center of mass (wrt center of chord)
  xp = ((xle_wbx*le_wbx + xte_wbx*te_wbx + x_ud*l_up+ x_ud*l_down)*t + xea*Ac)./(Ac+A_wbox);
  
  % dimensions of spar and wingbox
  % spar 
  bspar = Ac/2/hspar; % derived from the choice of bspar
  h = 0.14*chord;
  H = h+2*hspar;
  xspar = -0.15*chord;
  
  % wingbox (hollowed ellipsoid approx)
  a_ell = (0.75-0.1)*b;
  b_ell = h/2;
  a1 = a_ell-t;
  b1 = b_ell-t;
  
  % STIFFNESS PROPERTIES
  
  % AREA MOMENT OF INERTIAS (I_xx) (I-section + wingbox)
  % I-section
  Ispar = (tspar*h^3/12)+(bspar/12)*(H^3-h^3);
  % Ispar2 = 1/12*(bspar*H^3-h^3*(bspar-t)) % check
  % wingbox (hollowed ellipsoid approx)
  Iell = 2^4*pi/64*(a_ell*b_ell^3-a1*b1^3);
  % total beam section area moment of inertia [m^4] 
  I = Ispar+Iell;
  EI = E*I; 
  
  % TORSION CONSTANT of beam section [m^4] (only wingbox)
  U = pi*(a_ell+b_ell-t)*(1+0.258*(a_ell-b_ell)^2/(a_ell+b_ell-t)^2);
  A_ell = pi*(a_ell-0.5*t)*(b_ell-0.5*t); % mean of areas enclosed by outer
  % and inner boundaries, or (approximate) area within median boundary
  K = 4*A_ell^2*t/U;
  % J = pi*(a_ell^3*b_ell^3)/(a_ell^2+b_ell^2)*(1-(b1/b_ell)^4) 
  % there's a slight difference when using formula for constant thickness
  % hollowed ellipses (here K) and two ellipses one into another (here J)
  % doing J*G we have larger torsional stiffness so using K should be more
  % conservative (and also is a formula which is closer to reality)
  GK = G*K;
   
  % ROTARY INERTIAS (Jy = Iyy)
  
  % (spar)
  Ixx_spar = Ispar;
  Izz_spar = 1/12*(2*hspar*bspar.^3+h*tspar^3);
  Jspar = Ixx_spar+Izz_spar.^2;
  Jspar = rhop*Jspar+msparfun.*(xspar-xp).^2;
  
  % (ellipsoid)
  Ixx_ell = Iell;
  Izz_ell = 2^4*pi/64*(a_ell^3*b_ell-a1^3*b1); % [kgm]
  Jell = Ixx_ell+Izz_ell;
  Jell = rhop*Jell+m_wbox*(x_ud-xp).^2;
  
  % total  rotary inertia [kgm]
  Jp = Jell + Jspar; % [kgm]

  % aileron hinge axis mass properties
  md = 3*mhinge/l;            % hinge axis mass/unit span [kg/m]
  xd = 0.5*chord - chail;    % location of hinge axis axis [m]
  Jd = 0;                    % neglect rotary inertia [kgm]
  
  % equipment mass and inertia
  meq = 10;                 % equipment distributed mass [kg/m]
  xeq = 0;                  % equipment center of mass [m]
  Jeq = 3.5;                % equipment inertia [kgm]

  % center of distributed mass and rotary inertia
  xcm = (mp*xp+md*xd+meq*xeq) ./ (mp+md+meq);    % center of distributed mass
  Jcm = Jp+mp*(xp-xcm).^2 + Jd+md*(xd-xcm).^2 + Jeq+meq*(xeq-xcm).^2;   % rotary inertia wrt cm [kgm]

  % total distributed properties
  my = mp + md + meq;              % wing distributed mass/unit span [kg/m]
  s = xcm - xea;             % airfoil cm and ea separation [m]
  Jy = Jcm + my*s.^2;         % wing rotary inertia wrt ea/unit span [kgm]
  % Jy = mean(Jy);   % variations are not that big, it makes sense to assume 
  % this simplification

%  % stiffness properties
%   cplate = chord - chail;
%   I = (cplate*t^3/12); % beam section area moment of inertia [m^4]
%   K = (cplate*t^3/3);  % torsion constant of beam section [m^4]
%   EI = E*I;            % bending stiffness [Nm^2]
%   GK = G*K;            % torsional stiffness [Nm^2]
  
  % finally, setup wing with spanwise constant properties
  cb = ones(nsup, 1);
  geo = new_wing( y, xle*cb, xte*cb, ca, xea*cb, (xcm), ...
                  my*cb, (Jy), (EI), (GK)*cb, nelem );
 

  % add discrete mass definitions
  npm = size(dpm,1);
  for ki = 1:npm
      dJ = [dpm(ki,4) 0; 0 dpm(ki,5)];
    geo = attach_conmass(geo, dpm(ki,1), dpm(ki,2), dpm(ki,3), dJ);
  end

  % assemble FE matrices
  [M, K, Z] = assemble_structure(geo, B);
  M = Z' * M * Z;
  K = Z' * K * Z;

  % aileron aerodynamics
  [f, CRv, CRd] = assemble_aileron(geo);
  f = Z' * f;
  CRv = CRv * Z;

%   % construct interpolation tables for Qk
%   ndof = size(M, 1);
%   ktable = linspace(0.0, 2.0, 40);
%   nk = length(ktable);
%   Qtable = zeros(ndof,ndof,nk);
%   for ki = 1:nk
%     Qtable(:,:,ki) = Z' * assemble_aero(geo, ktable(ki)) * Z;
%   end
% 
%   Qip.ktab = ktable;
%   Qip.Qtab = Qtable;
%   Qip.bref = geo.bref;
%   Qip.rho = 1.225;

  % construct interpolation tables for Qk, also accounting for HTP
  ndof = size(M, 1);
  ktable = linspace(0.0, 2.0, 40);
  nk = length(ktable);
  Qtable = zeros(ndof,ndof,nk);
  
  % tail
  htp_c = 1.48;
  htp_l = 9.3/htp_c;
  xhtp = 8.5-b; % HTP center
  htp_b = 0.5*htp_c;
  T1node = [1 0 0; 0 0 0; -xhtp 0 1];
  T1 = [T1node zeros(3,3); zeros(3,3) T1node]; % not depending on red freq
  
  for ki = 1:nk
      Qac = assemble_aero(geo, ktable(ki)); % wing
      
      khtp = ktable(ki)*htp_b/geo.bref; % scaling reduced freq for tail
      Qhtp = beam_amatrix(khtp, htp_l, htp_b, 0);
      T2node = [1 0 -xhtp; 0 0 0; 0 0 (1+1i*ktable(ki)*xhtp/b)];
      % last elem of T2node also account for delay related to angle of
      % attack between main wing and tail plane 
      % T2node = [1 0 -xhtp; 0 0 0; 0 0 1]; % not accounting for delay
      T2 = [T2node, zeros(3,3); zeros(3,3), T2node];
      
      % r position (where adding contribution of tail to the complete
      % aircraft) should be in the middle of the wing. 
      % r is an interval of 6 dofs
      n = (1+nelem)*3/2; % choose odd nelem
      r = (n-2):(n+3);
      Qac(r,r) = Qac(r,r) + T1*Qhtp*T2;
      Qtable(:,:,ki) = Z'*Qac*Z;
      %Qtable(:,:,ki) = Z' * assemble_aero(geo, ktable(ki)) * Z;
  end

  Qip.ktab = ktable;
  Qip.Qtab = Qtable;
  Qip.bref = geo.bref;
  Qip.rho = 1.225;

end

function [geo] = new_wing(y, xle, xte, ca, xea, xcm, my, Jy, EI, GK, nelem)
% [geo] = new_wing(y, xle, xte, ca, xea, xcm, my, Jy, EI, GK, nelem)
%
% Defines geometric and structural properties for a wing model
% based on beam theory. This is the most general version which
% assumes that all properties vary along the span.
%
% Input arguments:
% y     : Vector of spanwise positions at which the remaining
%         properties are given. These do not necessarily coincide
%         with the beam FEM nodes.
% xle   : x-coordinates of the leading edge at y
% xte   : x-coordinates of the trailing edge at y
% ca    : aileron chord at y (ca(i) = 0 if none present at y(i))
% xea   : x-coordinate of the elastic axis at y
% xcm   : x-coordinate of the center of mass at y
% my    : mass per length at y [kg/m]
% Jy    : torsional inertia per length at y [kgm?/m]
% nelem : number of finite elements required
%
% Output struct:
% geo.y    : spanwise position of the finite element nodes (vector)
% geo.x    : chordwise position of the midchord points (vector)
% geo.b    : semichord at FE nodes (vector)
% geo.bref : reference semichord (scalar)
% geo.ba   : aileron semichord at FE nodes, zero if no aileron (vector)
% geo.xea  : chordwise position of elastic axis (vector)
% geo.xcm  : chordwise position of center of mass (vector)
% geo.my   : mass per length at nodes [kg/m] (vector)
% geo.jy   : torsional inertia per length at nodes [kgm] (vector)
% geo.ei   : bending stiffness at nodes [Nm^2] (vector)
% geo.gk   : torsional stiffness at nodes [Nm^2] (vector)
% geo.rhof : fluid density, default 1.225 [kg/m^3] (scalar)
%
% (c) 2004-2016 Dan Borglund <dodde@kth.se> and David Eller <dlr@kth.se>

  % first, compute planform values for spanwise node
  % positions and chord dimensions
  yt = linspace(min(y), max(y), nelem+1);
  xl = interp1(y, xle, yt, 'spline');
  xt = interp1(y, xte, yt, 'spline');
  % vb = 0.5*(xt - xl);

  % fill struct with planform definition
  geo.y = yt;
  geo.x = 0.5*(xl + xt);
  geo.b = 0.5*(xt - xl);

  % print area (checking)
  area = trapz(yt, xt-xl);
  fprintf(1, 'Wing area: %.3f m^2\n', area);

  % semichord used as reference value
  geo.bref = 0.5*area / (max(yt) - min(yt));

  % compute aileron semichord vector (no interpolation)
  geo.ba = 0.5 * interp1(y, ca, yt, 'nearest');

  % interpolate mass properties
  geo.xea = interp1(y, xea, yt, 'spline');
  geo.xcm = interp1(y, xcm, yt, 'spline');
  geo.my = interp1(y, my, yt, 'spline');
  geo.jy = interp1(y, Jy, yt, 'spline');

  % the default wing does not have any lumped masses
  geo.cmi = [];
  geo.cdm = [];
  geo.cmx = [];
  geo.cmy = [];
  geo.cdj = {};

  % print total mass and inertia (checking)
  mass = trapz(yt, geo.my);
  tinert = trapz(yt, geo.jy);
  fprintf(1, 'Wing mass: %.3f kg\n', mass);
  fprintf(1, 'Torsional inertia: %.3f kgm^2\n', tinert);

  % interpolate structural properties
  geo.ei = interp1(y, EI, yt, 'spline');
  geo.gk = interp1(y, GK, yt, 'spline');
end

function [geo] = attach_conmass(geo, dm, x, y, dJ)
% [geo] = attach_conmass(geo, dm, x, y, dJ)
%
% Update the wing struct geo with a concentrated mass
% specification.
%
% geo:  the wing struct, modified copy is returned
% dm:   concentrated mass, scalar [kg]
% x:    longitudinal position of the concentrated mass [m]
% y:    spanwise position of the concentrated mass [m]
% dJ:   2x2 rotational inertia matrix [kgm²] (optional argument)
%
% (c) 2004-2016 Dan Borglund <dodde@kth.se> and David Eller <dlr@kth.se>

  % locate attachment node
  vklo = find(geo.y <= y);
  vkhi = find(geo.y >= y);
  klo = vklo(end);
  khi = vkhi(1);

  % pick the nearest node
  dlo = (geo.x(klo) - x)^2 + (geo.y(klo) - y)^2;
  dhi = (geo.x(khi) - x)^2 + (geo.y(khi) - y)^2;
  if dlo < dhi
    inode = klo;
  else
    inode = khi;
  end

  % update wing definition
  geo.cmi = [geo.cmi; inode];
  geo.cdm = [geo.cdm; dm];
  geo.cmx = [geo.cmx; x];
  geo.cmy = [geo.cmy; y];

  if nargin > 4
    geo.cdj = cat(1, geo.cdj, {dJ});
  else
    dJ = zeros(2,2);
    geo.cdj = cat(1, geo.cdj, {dJ});
  end

end

function [M, K, Z] = assemble_structure(geo, B)
% [M, K, Z] = assemble_structure(geo, B)
%
% Assemble global mass and stiffness matrices for the beam FE
% model. The input argument geo is a struct witch defines the
% wing geometry and properties.
%
% geo :  description of structure (struct)
% dsy :  aeroelastic system (struct)
% B :    optional argument, linear constraints on DOFS
%
% M  :  unconstrained mass matrix
% K  :  unconstrained stiffness matrix
% Z  :  subspace projection for constraints or modal subspace
%
% (c) 2004-2016 Dan Borglund <dodde@kth.se> and David Eller <dlr@kth.se>

  % calculate dimensions
  nelem = length(geo.y) - 1;
  ndof = 3 * (nelem + 1);
  M = zeros(ndof, ndof);
  K = zeros(ndof, ndof);

  for k = 1:nelem

    % average element values from nodal properties
    ly = geo.y(k+1) - geo.y(k);
    lx = geo.x(k+1) - geo.x(k);
    le = sqrt(lx^2 + ly^2);
    xcm = 0.5 * ( geo.xcm(k+1) + geo.xcm(k));
    xea = 0.5 * ( geo.xea(k+1) + geo.xea(k));
    se = xcm - xea;
    mye = 0.5 * ( geo.my(k+1) + geo.my(k) );
    Jye = 0.5 * ( geo.jy(k+1) + geo.jy(k) );
    EIe = 0.5 * ( geo.ei(k+1) + geo.ei(k) );
    GKe = 0.5 * ( geo.gk(k+1) + geo.gk(k) );

    % fetch element matrices
    Ke = beam_kmatrix(le, EIe, GKe);
    Me = beam_mmatrix(le, se, mye, Jye);

    % find index range where the element matrices are
    % assembled into the global system matrices
    first = 3*(k-1) + 1;
    last  = first + 5;
    rng = first:last;

    % assemble
    M(rng,rng) = M(rng,rng) + Me;
    K(rng,rng) = K(rng,rng) + Ke;

  end

  % update mass matrix with concentrated masses
  nmc = length(geo.cdm);
  for k = 1:nmc

    % extract specification from struct
    inode = fix(geo.cmi(k));
    dm = geo.cdm(k);
    dx = geo.cmx(k) - geo.xea(inode);
    dy = geo.cmy(k) - geo.y(inode);
    dJ = geo.cdj{k};

    % rotational inertia relative to inode
    Jn = zeros(2,2);
    Jn(1,1) = dJ(1,1) + dm*dy^2;
    Jn(1,2) = dJ(1,2) + dm*dx*dy;
    Jn(2,1) = dJ(2,1) + dm*dx*dy;
    Jn(2,2) = dJ(2,2) + dm*dx^2;

    % nodal mass matrix
    Mn = zeros(3,3);

    % mass and static moments
    Mn(1,1) = dm;
    Mn(1,2) = -dm*dy;
    Mn(1,3) = -dm*dx;
    Mn(2,1) = Mn(1,2);
    Mn(3,1) = Mn(1,3);

    % rotational inertias
    Mn(2:3,2:3) = Jn;

    % find DOF indices corresponding to inode
    first = 3*(inode - 1) + 1;
    last = first + 2;
    rng = first:last;

    % check input
    if first < 1 || last > ndof
      error('Concentrated mass node index out of range.');
    end

    % update global mass matrix
    M(rng,rng) = M(rng,rng) + Mn;

  end  % concentrated masses

  % set projection
  if nargin == 2 && ~isempty(B)
    Z = null(B);
  else
    Z = eye(ndof);
  end

end

function [Ak] = assemble_aero(geo, k)
% [Ak] = assemble_aero(geo, k)
%
% Assemble global unsteady aerodynamic load matrix Ak
% for the wing defined in geo and reduced frequency k.
%
% (c) 2004-2016 Dan Borglund <dodde@kth.se> and David Eller <dlr@kth.se>

  % calculate dimensions
  nelem = length(geo.y) - 1;
  ndof = 3 * (nelem + 1);
  Ak = zeros(ndof, ndof);

  % extract node arrays
  vy = geo.y;
  vb = geo.b;
  vxea = geo.xea;
  vx = geo.x;

  for ki = 1:nelem

    % average element values from neighbor nodal properties
    le = vy(ki+1) - vy(ki);
    b = 0.5*(vb(ki+1) + vb(ki));
    xea = 0.5*(vxea(ki+1) + vxea(ki));
    xmc = 0.5*(vx(ki+1) + vx(ki));
    xa = xea - xmc;

    % reduced frequency for this element
    rfe = k * b / geo.bref;

    % fetch and scale element matrix
    Ae = beam_amatrix(rfe, le, b, xa);

    % find index range where the element matrices are
    % assembled into the global system matrices
    first = 3*(ki-1) + 1;
    last  = first + 5;
    rng = first:last;

    % assemble
    Ak(rng,rng) = Ak(rng,rng) + Ae;
  end

end

function [f,CRv,CRd] = assemble_aileron(geo)
% [f,CRv,CRd] = assemble_aileron(geo)
%
% Assemble aileron load vectors. With the results, the rolling moment R
% caused by a deformation vector v and an aileron deflection delta
% (in radian) can be computed as:
%
%        R = q*( Rv*v + Rd*delta )
%
% where v is assumed to be a column vector.
%
% geo:  struct describing discretization and aileron geometry
% f:    q*f*delta gives the nodal loads as a function of aileron
%       deflection delta (column vector)
% CRv:  rolling moment coefficient due to deformation (row vector)
% CRd:  rolling moment due to aileron deflection (scalar)
%
% (c) 2004-2016 Dan Borglund <dodde@kth.se> and David Eller <dlr@kth.se>

  % extract section dimensions
  l = max(geo.y) - min(geo.y);    % span
  Sref = 2*trapz(geo.y, geo.b);   % wing area

  % calculate sizes, reserve space
  nelem = length(geo.y) - 1;
  ndof = 3 * (nelem + 1);
  Rv = zeros(1, ndof);
  f = zeros(ndof, 1);
  Rd = 0;

  % assemble from root to tip
  for ki=1:nelem

    b = geo.b(ki);                        % local semichord
    a = geo.xea(ki) / geo.bref;           % nondimensional location of cr
    xd = b - 2*geo.ba(ki);
    e = xd / geo.bref;                    % nondimensional location of aileron hinge
    y1e = geo.y(ki);                      % span location of first node
    elene = geo.y(ki+1) - y1e;            % element width
    elene2 = elene^2;

    % static aerodynamic coefficients
    CLa = 2*pi;
    CLd = 2*(acos(e)+sqrt(1-e^2));
    % CMa = pi*(1/2+a);
    CMd = (CLd*a-(e*sqrt(1-e^2)-acos(e)))/2;

    % compute 6x1 element vectors
    fe = zeros(6,1);
    Rve = zeros(1,6);

    fe(1) = (1/2)*(2*b)*CLd*elene;
    fe(2) = (1/12)*(2*b)*CLd*elene2;
    fe(3) = (1/2)*(2*b)^2*CMd*elene;
    fe(4) = (1/2)*(2*b)*CLd*elene;
    fe(5) = -(1/12)*(2*b)*CLd*elene2;
    fe(6) = (1/2)*(2*b)^2*CMd*elene;

    Rve(1) = 0;
    Rve(2) = 0;
    Rve(3) = (2*b)*CLa*(y1e*(1/2)*elene+(1/6)*elene2);
    Rve(4) = 0;
    Rve(5) = 0;
    Rve(6) = (2*b)*CLa*(y1e*(1/2)*elene+(1/3)*elene2);

    Rde = (2*b)*CLd*elene*(y1e+(1/2)*elene);

    % add element vector Rve to global vector Rv
    first = 3*(ki-1) + 1;
    last  = first + 5;
    rng = first:last;
    Rv(rng) = Rv(rng) + Rve;

    % add element vector fe and coefficient Rde if aileron element
    if geo.ba(ki) > 0
      f(rng) = f(rng) + fe;
      Rd = Rd + Rde;
    end
  end

  % convert to coefficients
  CRd = Rd / (Sref*l);
  CRv = Rv / (Sref*l);

end

function [Me] = beam_mmatrix(le, se, mye, Jye)
% [Me] = beam_mmatrix(le, se, mye, Jye)
%
% Construct the element mass matrix for a beam element with element
% length le, mass per length mye, and torsional intertia per length Jye.
% se is the distance between center of mass and elastic axis xcm - xea.
%
% (c) 2004-2016 Dan Borglund <dodde@kth.se> and David Eller <dlr@kth.se>

  Me = zeros(6,6);
  le2 = le^2;
  le3 = le^3;

  Me(1,1) = (13/35)*mye*le;
  Me(1,2) = (11/210)*mye*le2;
  Me(1,3) = -(7/20)*mye*se*le;
  Me(1,4) = (9/70)*mye*le;
  Me(1,5) = -(13/420)*mye*le2;
  Me(1,6) = -(3/20)*mye*se*le;
  Me(2,2) = (1/105)*mye*le3;
  Me(2,3) = -(1/20)*mye*se*le2;
  Me(2,4) = (13/420)*mye*le2;
  Me(2,5) = -(1/140)*mye*le3;
  Me(2,6) = -(1/30)*mye*se*le2;
  Me(3,3) = (1/3)*Jye*le;
  Me(3,4) = -(3/20)*mye*se*le;
  Me(3,5) = (1/30)*mye*se*le2;
  Me(3,6) = (1/6)*Jye*le;
  Me(4,4) = (13/35)*mye*le;
  Me(4,5) = -(11/210)*mye*le2;
  Me(4,6) = -(7/20)*mye*se*le;
  Me(5,5) = (1/105)*mye*le3;
  Me(5,6) = (1/20)*mye*se*le2;
  Me(6,6) = (1/3)*Jye*le;

  Me = Me + triu(Me,1)'; % exploit symmetry for the lower triangular part
end

function [Ae] = beam_amatrix(k, le, b, xa)
% [Ae] = beam_amatrix(k, le, b, xa)
%
% Construct unsteady aerodynamic load matrix Ae for a
% beam element with length le and semichord b. xa is the
% local position of the elastic axis relative the the
% midchord (xea - xmid).
%
% (c) 2004-2016 Dan Borglund <dodde@kth.se> and David Eller <dlr@kth.se>

  % nondimensional distance xea - midchord
  a = xa/b;

  % Theodorsen's function
  if k == 0
    Ck = 1;
  else
    H0 = besselh(0,2,k);
    H1 = besselh(1,2,k);
    Ck = H1/(H1+i*H0);
  end

  % Theodorsen's frequency-domain aerodynamic coefficients
  Lh = -pi*((i*k)^2+2*Ck*(i*k));
  La = -pi*(a*(i*k)^2-(1+2*(1/2-a)*Ck)*(i*k)-2*Ck);
  Mh = -(pi/2)*(a*(i*k)^2+2*(1/2+a)*Ck*(i*k));
  Ma = -(pi/2)*((1/8+a^2)*(i*k)^2+((1/2-a)-2*(1/2+a)*(1/2-a)*Ck)*(i*k) -2*(1/2+a)*Ck);

  Ae = zeros(6,6);
  le2 = le^2;
  le3 = le^3;

  Ae(1,1) = (13/35)*(2*b)*(Lh/b)*le;
  Ae(1,2) = (11/210)*(2*b)*(Lh/b)*le2;
  Ae(1,3) = (7/20)*(2*b)*La*le;
  Ae(1,4) = (9/70)*(2*b)*(Lh/b)*le;
  Ae(1,5) = -(13/420)*(2*b)*(Lh/b)*le2;
  Ae(1,6) = (3/20)*(2*b)*La*le;

  Ae(2,1) = (11/210)*(2*b)*(Lh/b)*le2;
  Ae(2,2) = (1/105)*(2*b)*(Lh/b)*le3;
  Ae(2,3) = (1/20)*(2*b)*La*le2;
  Ae(2,4) = (13/420)*(2*b)*(Lh/b)*le2;
  Ae(2,5) = -(1/140)*(2*b)*(Lh/b)*le3;
  Ae(2,6) = (1/30)*(2*b)*La*le2;

  Ae(3,1) = (7/20)*(2*b)^2*(Mh/b)*le;
  Ae(3,2) = (1/20)*(2*b)^2*(Mh/b)*le2;
  Ae(3,3) = (1/3)*(2*b)^2*Ma*le;
  Ae(3,4) = (3/20)*(2*b)^2*(Mh/b)*le;
  Ae(3,5) = -(1/30)*(2*b)^2*(Mh/b)*le2;
  Ae(3,6) = (1/6)*(2*b)^2*Ma*le;

  Ae(4,1) = (9/70)*(2*b)*(Lh/b)*le;
  Ae(4,2) = (13/420)*(2*b)*(Lh/b)*le2;
  Ae(4,3) = (3/20)*(2*b)*La*le;
  Ae(4,4) = (13/35)*(2*b)*(Lh/b)*le;
  Ae(4,5) = -(11/210)*(2*b)*(Lh/b)*le2;
  Ae(4,6) = (7/20)*(2*b)*La*le;

  Ae(5,1) = -(13/420)*(2*b)*(Lh/b)*le2;
  Ae(5,2) = -(1/140)*(2*b)*(Lh/b)*le3;
  Ae(5,3) = -(1/30)*(2*b)*La*le2;
  Ae(5,4) = -(11/210)*(2*b)*(Lh/b)*le2;
  Ae(5,5) = (1/105)*(2*b)*(Lh/b)*le3;
  Ae(5,6) = -(1/20)*(2*b)*La*le2;

  Ae(6,1) = (3/20)*(2*b)^2*(Mh/b)*le;
  Ae(6,2) = (1/30)*(2*b)^2*(Mh/b)*le2;
  Ae(6,3) = (1/6)*(2*b)^2*Ma*le;
  Ae(6,4) = (7/20)*(2*b)^2*(Mh/b)*le;
  Ae(6,5) = -(1/20)*(2*b)^2*(Mh/b)*le2;
  Ae(6,6) = (1/3)*(2*b)^2*Ma*le;
end

function [fg] = assemble_gustcolumn(geo, k)
% [fg] = assemble_gustcolumn(geo, k)
%
% Assemble element gust load vectors into a global gust column vector.
% rho*uoo*wg * fg(k) = 2*qoo * fg(k) * wg/uoo yields the vector of dimensional
% forces and moment which is the right-hand side for a harmonic vertical gust
% with vertical speed wg. Here, the reduced frequency is determined from the
% gust wavelength Lw by k = 2*pi*b / Lw
%
% (c) 2016 David Eller <dlr@kth.se>

  % calculate dimensions
  nelem = length(geo.y) - 1;
  ndof = 3 * (nelem + 1);
  fg = zeros(ndof, 1);

  % extract node arrays
  vy = geo.y;
  vb = geo.b;
  vxea = geo.xea;
  vx = geo.x;

  for ki = 1:nelem

    % average element values from neighbor nodal properties
    le = vy(ki+1) - vy(ki);
    b = 0.5*(vb(ki+1) + vb(ki));
    xea = 0.5*(vxea(ki+1) + vxea(ki));
    xmc = 0.5*(vx(ki+1) + vx(ki));
    xa = xea - xmc;

    % reduced frequency for this element
    rfe = k * b / geo.bref;

    % fetch element vector
    fge = beam_gustcolumn(rfe, le, b, xa);

    % assemble
    first = 3*(ki-1) + 1;
    last  = first + 5;
    rng = first:last;
    fg(rng) = fg(rng) + fge;
  end

end

function [fge] = beam_gustcolumn(k, le, b, xa)
% [fge] = beam_gustcolumn(k, le, b, xa)
%
% Construct harmonic gust load vector fg for a beam element with length le and
% semichord b. xa is the local position of the elastic axis relative the the
% midchord (xea - xmid). In order to obtain dimensional forces and moments,
% multiply fg by rho*uoo*wg, where rho is the air density, uoo the speed of
% flight and wg is vertical gust velocity.
%
% (c) 2016 David Eller <dlr@kth.se>

  % nondimensional distance xea - midchord
  a = xa/b;

  % evaluate for positive frquency first
  ka = abs(k);

  % Theodorsen's function
  if ka == 0
    Ck = 1;
  else
    H0 = besselh(0,2,ka);
    H1 = besselh(1,2,ka);
    Ck = H1 / (H1 + 1i*H0);
  end

  % Sear's function
  J0 = besselj(0, ka);
  J1 = besselj(1, ka);
  Sk = Ck*(J0 - 1i*J1) + 1i*J1;

  % make sure this works when called with negative frequency
  if k < 0
    Sk = conj(Sk);
  end

  % Sear's frequency-domain aerodynamic coefficients
  Lg = 2*pi*b*Sk;
  Mg = b*(0.5 + a)*Lg;

  fge = zeros(6,1);
  le2 = le^2;

  fge(1) = (1/2)*Lg*le;
  fge(2) = (1/12)*Lg*le2;
  fge(3) = (1/2)*Mg*le;

  fge(4) = (1/2)*Lg*le;
  fge(5) = -(1/12)*Lg*le2;
  fge(6) = (1/2)*Mg*le;

end

function [Ke] = beam_kmatrix(le, EIe, GKe)
% [Ke] = beam_kmatrix(le, EIe, GKe)
%
% Construct beam element stiffness matrix for an element
% with length le and stiffness properties EIe, GKe.
%
% (c) 2004-2016 Dan Borglund <dodde@kth.se> and David Eller <dlr@kth.se>

  Ke = zeros(6,6);
  le2 = le^2;
  le3 = le^3;

  Ke(1,1) = 12*EIe/le3;
  Ke(1,2) = 6*EIe/le2;
  Ke(1,3) = 0;
  Ke(1,4) = -12*EIe/le3;
  Ke(1,5) = 6*EIe/le2;
  Ke(1,6) = 0;
  Ke(2,2) = 4*EIe/le;
  Ke(2,3) = 0;
  Ke(2,4) = -6*EIe/le2;
  Ke(2,5) = 2*EIe/le;
  Ke(2,6) = 0;
  Ke(3,3) = GKe/le;
  Ke(3,4) = 0;
  Ke(3,5) = 0;
  Ke(3,6) = -GKe/le;
  Ke(4,4) = 12*EIe/le3;
  Ke(4,5) = -6*EIe/le2;
  Ke(4,6) = 0;
  Ke(5,5) = 4*EIe/le;
  Ke(5,6) = 0;
  Ke(6,6) = GKe/le;

  Ke = Ke + triu(Ke,1)'; % exploit symmetry for the lower triangular part
end
