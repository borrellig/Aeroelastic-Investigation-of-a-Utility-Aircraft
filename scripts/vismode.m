function vismode(v,np)
% vismode(geo,v,np)
%
% Visualize a single modeshape using a matlab 'movie'. Note that the
% deformation vector v must be in the original space used to define
% the geometry, not a projected deformation.
%
% v :     modeshape to visualize (vector)
% np :    number of cycles to plot (optional argument)
%
% (c) 2016 Dan Borglund <dodde@kth.se> and David Eller <dlr@kth.se>

  % constant dimensions
  b = 0.15;
  l = 1.2;
  xea = -0.03;

  % generate np periods of oscillation with ns snapshots/period
  % if the display is too fast, use more snapshots/period, i.e.
  % increase ns to 200 or more
  ns = 800;

  % show 5 periods by default, if no other arguments are given
  if nargin < 3
    np = 5;
  end

  % check dimension of v
  ndof = length(v);
  if numel(v) ~= ndof
    error('vismode: Can only plot one modeshape at a time.');
  elseif mod(ndof,3) ~= 0
    error('vismode: Number of DOFs is not a multiple of 3.');
  end

  % fetch dimensions and spanwise coordinates
  nnode = fix(ndof/3);
  nshape = [1, nnode];
  yp = linspace(0.0, l, nnode);

  % Extract bending and torsion dofs section by section
  w = zeros(nshape);
  t = zeros(nshape);
  for k = 1:nnode
    w(k) = v(1+3*(k-1)); % nodal deflection
    t(k) = v(3+3*(k-1)); % nodal twist
  end

  % compute le and te coordinates at spanwise locations of nodes
  xte = b * ones(nshape);
  xle = -xte;
  yte = yp;
  yle = yp;
  zle = w + (b + xea)*ones(nshape) .* t;
  zte = w - (b - xea)*ones(nshape) .* t;

  % visualize using the mesh feature
  X = [xle;xte];
  Y = [yle;yte];

  % scale maximum amplitude to b/2
  maxdef = max( abs([zle zte]));
  Z = (b/2)*[zle;zte] / maxdef;

  % precompute phase angles
  nframes = np*ns;
  vphi = linspace(0.0, np*2*pi, nframes);

  % render first frame
  kj = 1;
  clf;
  hold on;
  phi = vphi(kj) + angle(Z);
  h = mesh(X, Y, abs(Z) .* sin(phi));
  plot3([0 1.2*b],[0 0],[0 0],'b');
  text(1.4*b,0,0,'x');
  plot3([0 0],[0 0],[0 0.5*b],'b');
  text(0,0,0.7*b,'z');
  axis([-2*b 2*b 0 l -2*b 2*b]);
  axis off;
  view (130,20);
  hold off;
  drawnow;

  % render remaining frames by replacing deflection data
  % in plot handle only; do not redraw everything
  for kj = 2:nframes
    phi = vphi(kj) + angle(Z);
    set(h, 'ZData', abs(Z) .* sin(phi));
    drawnow;
  end

end
