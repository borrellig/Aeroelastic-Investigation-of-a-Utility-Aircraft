function [y,we,thetae,Mb,Mt,S_normal,S_shear] = stresses(v,EI,GK,l,nelem,H,I,K)

% nnodes = nelem+1;
l_e = l/nelem;
% basis function (y accounts for eps, see p.27)
Nw = @(y) [(1-3*y.^2+2*y.^3) 
    l_e*(y-2*y.^2+y.^3) 
    3*y.^2-2*y.^3 
    l_e*(-y.^2+y.^3)];
Ntheta = @(y) [1-y; y];

% derivatives related to the change of coordinate (y-to-eps)
depsdy = 1/l_e;
ddepsdy = 1/l_e^2;

% evaluation of basis functions at the nodes
% y = yn(1:end-1);
y = l/2*linspace(-1,1,nelem);
Nwy = Nw(y);
Nthetay = Ntheta(y);
n = length(y);

% interpolation within the element (displacement and twist in element)
for i = 1:length(y)-1
    we(i) = sum(Nwy(1:2,i)'*v((3*(i-1)+1):(3*(i-1)+2))+Nwy(3:4,i)'*v((3*i+1):(3*i+2)));
    thetae(i) = Nthetay(1,i)*v(3*(i-1)+3)+Nthetay(2,i)*v(3*i+3); % radians
end

%% bending and torsion moments using finite differences on the basis
h = 0.01;

for i = 1:length(y)-1
    % w''
    Nw2(:,i) = (Nw(y(i)+h)-2*Nw(y(i))+Nw(y(i)-h))/h^2;
    w2(i) = sum(Nw2(1:2,i)'*v((3*(i-1)+1):(3*(i-1)+2))+Nw2(3:4,i)'*v((3*i+1):(3*i+2)));
    % theta'
    Ntheta1(:,i) = 0.5*(Ntheta(y(i)+h)-Ntheta(y(i)-h))/h;
    theta1(i) = Ntheta1(1,i)*v(3*(i-1)+3)+Ntheta1(2,i)*v(3*i+3);
    
    % bending moment
    if length(EI) == 1
        Mb(i) = EI*w2(i)*ddepsdy;
    else
        Mb(i) = EI(i)*w2(i)*ddepsdy;
    end
    % twist moment
    Mt(i) = GK*theta1(i)*depsdy;
    % remember to multiply by the derivative of epsilon wrt y, i.e.
    % eps = (y-y1)/l_e, d(eps)/dy = 1/l_e
end

%% stresses (normal and shear)
% first defining some locations (arbitrary) along the thickness where to
% analyze the stresses
nz = 5;
z = H/2*linspace(-1,1,nz);

% stresses evaluation at the chosen locations
for i = 1:length(y)-1
    % normal stress distribution
    if length(EI) == 1
        S_normal(:,i) = -Mb(i)/I*z; % negative sign because in tension stresses
        % are positive in the lower part of the wing
    else
        S_normal(:,i) = -Mb(i)/I(i)*z; %N/m^2
    end
    
    % shear stress distribution
    S_shear(:,i) = -2*Mt(i)/K*z; % thin-walled open section
end
   
y = y-l_e; % in order to position center of fullspan

%% specification for safety factors (Al2024-T3)
Sn_ult = 483E6; % MPa, ultimate tensile strength
Ss_ult = 283E6; % MPa, shear strength

Sn_max = abs(max(max(S_normal)));
Ss_max = abs(max(max(S_shear)));

fprintf('\nMax. Normal stress = %.2f MPa\n',Sn_max/1E6)
fprintf('\nMax. Shear stress = %.2f MPa\n',Ss_max/1E6)

%safety coefficients (need to be above 1.5)
Snratio = Sn_ult/Sn_max;
Ssratio = Ss_ult/Ss_max;

fprintf('\nMin. Safety coefficient (normal) = %.2f\n',Snratio)
fprintf('Min. Safety coefficient (shear) = %.2f\n',Ssratio)


