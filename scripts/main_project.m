clear all
close all

fprintf('<strong> ------------------------------------------- \n</strong>')
format shortG
% certification limits
fprintf('<strong> \nCertification limits: \n</strong>')
v_ne = 200*0.51; % m/s
v_cruise = 150*0.51;
mass_max = 5670; % kg
fuel_capacity = 1179; % kg

fprintf('V_ne = %.2f m/s \n',v_ne)
fprintf('Max take-off weight: %.2f kg\n',mass_max)

%% setup geometry and structural properties 
fprintf('<strong> \nGeometry and structural properties: \n</strong>')

nelem = 99; % choose odd nelem (related to the introduction of the HTP)
nnodes = nelem+1;
ndof = 3*nnodes;

l = 19.8;
S = 39;
b = 0.5*S/l;
ba = (1-0.8)*b;
mhinge = 0;
rhop = 2780;
E = 73.1E9; % 10.6E3 ksi (1 ksi = 6894757.2931783 Pa)
G = 28E9; % 4.06E3 ksi

% Dimensions for the most critical case (n = 4.4)
yn = l/2*linspace(-1,1,nnodes);
Ac = 0.01;

% constant spar
% ts = 0.001; 
% Acfun = @(y) Ac-0*abs(yn); 

% tapered spar
ts = 0.0009; % tapered spar
Acfun = @(y) Ac-0.001*abs(yn);

Ac = Acfun(yn);
hspar = 0.01;
bspar = Ac/2/hspar;

% Tail specifications
htp_c = 1.48;
htp_S = 9.3;
htp_l = htp_S/htp_c;
xhtp = 8.5-b; % HTP center
htp_b = 0.5*htp_c;

% concentrated masses
% dpm = [];
dpm = createdpm(b);

% constraints + retrieve system matrices (FULLSPAN)
B = zeros(3,ndof);
B(1,1:3:end) = 1;
B(2,2:3:end) = 1;
B(2,1:3:end) = yn;
B(3,3:3:end) = 1;

% clamped wing (!!! tests)
% B = zeros(3,ndof);
% yn = linspace(0,l,nnodes);
% B(1,1) = 1;
% B(2,2) = 1;
% B(3,3) = 1;

% Z = null(B);
% K = Z'*K*Z;
% M = Z'*M*Z;

% check of total mass after adding concentrated masses
vacc=zeros(3*nnodes,1);
vacc(1:3:ndof)=1;
% [M,~,~,~,~,~,~,~,~] = otterwing([], l, b, ts, ba, mhinge, rhop, E, G, nelem, dpm, Ac, hspar);
M1 = otterwing([], l, b, ts, ba, mhinge, rhop, E, G, nelem, dpm, Ac,hspar);
% total_mass = vacc'*M1*vacc;
total_mass = B(1,:)*M1*B(1,:)'; % for fullspan wing

% vibration eigenvalue problem
% [M,K,Z,Qip,f,CRv,CRd,EI,GK] = otterwing(B,l,b,ts,ba,mhinge,rhop,E,G,nelem,dpm,Ac, hspar);
[M,K,Z,Qip,f,CRv,CRd,EI,GK] = otterwing(B,l,b,ts,ba,mhinge,rhop,E,G,nelem,dpm,Ac, hspar);
% check if rigid body motion
% norm(Z*K*Z'*B') % should be equal to zero for rigid body motion
[V,D] = eig(K,M); % visualize with vismode
[ev,ind] = sort(diag(D));
V = V(:,ind);
check_freq = sqrt(ev(1:8))/2/pi;

fprintf('\nTotal weight: %.2f kg\n',total_mass)
fprintf('Total payload allowed: %.2f kg\n',mass_max-total_mass)
if total_mass>mass_max
    fprintf('Error: weight above limitations')
end

%% TASK3)
fprintf('<strong> \nTASK3) Aeroelastic Trim: \n</strong>')
% static aeroelastic trim in order to find (alpha,d_e,v) for different
% flight loads (including inertial relief) + plot spanwise distribution of
% shear and normal stresses at most critical flight conditions
% 
g = 9.81;
n = 4.4;
% n = -1.8;
% n = 1;
CLmax = 1.5;

%%% IMPORTANT, if wanting to use CL_de and Cm_de from tables it's important
%%% to use the adyn matrix Qip.Qtable without wing contribution, only in
%%% that way it is possible to obtain values of CLalpha and Cmalpha close
%%% to the one given

if n ==  1
    u = v_cruise; % cruise velocity used for trim (at n = 1)
else
    % increase speed when higer load factors (reasonable to reach them)
    u = v_ne;
end
q = 0.5*Qip.rho*u^2;

% CLalpha = 5.66;
CL_de = 0.608; % positive deflection (dwns) generating lift up
% Cm_alpha = -1.31;
Cm_de = -1.74; % positive deflection generating pitch down (e.g. descending)

[M,K,Z,Qip,f,CRv,CRd,EI,GK] = otterwingTrim([],l,b,ts,ba,mhinge,rhop,E,G,nelem,dpm,Ac,hspar);
Q0=ipolQk(Qip,0.0);
e1 = zeros(ndof,1);
e3 = e1;
e1(1:3:end) = 1;
e3(3:3:end) = 1;
Z = null(B);

% Ael. feedback
lhs = [Z'*(K-q*Q0)*Z  -Z'*q*Q0*e3 zeros(size(Z,2),1)
    e1'*q*Q0*Z  e1'*q*Q0*e3  q*htp_S*CL_de
    e3'*q*Q0*Z  -q*e3'*Q0*e3  q*htp_S*htp_c*Cm_de];
rhs = [-Z'*M*e1*n*g
    e1'*M*e1*n*g
    0];

x = lhs\rhs;

v = Z*x(1:end-2)+e3*x(end-1);
alpha = x(end-1);
d_e = x(end);
Lift = q*(Q0*v)'*e1+q*htp_S*CL_de*d_e;

fprintf('\nAlpha_trim = %.2f deg \n',rad2deg(alpha))
fprintf('Elevator deflection = %.2f deg \n',rad2deg(d_e))
fprintf('CLtrim = %.2f (needs to be < 1.5)\n',abs(e1'*Q0*e3/S*alpha))
fprintf('CL_elevator = %.2f (needs to be < 1.5)\n',abs(CL_de*d_e))

figure
subplot(2,1,1)
title('Vertical displacement')
xlabel('y [m]')
ylabel('w [m]')
hold on
plot(yn,v(1:3:end),'o-')
grid on

subplot(2,1,2)
title('Torsion angle')
xlabel('y [m]')
ylabel('theta [rad]')
hold on
plot(yn,v(3:3:end),'o-')
grid on

%% TASK2) 
fprintf('<strong> \nTASK2) Stress analysis: \n</strong>')
% stresses analysis
% 
fprintf('\nLoad factor = %.1f, velocity = %.2f m/s\n',n,u)
h = 0.14*2*b;
H = h+2*hspar;
v = v(1:end-1);
[y,we,thetae,Mb,Mt,S_normal,S_shear] = stresses(v,EI,GK,l,nelem,H,EI/E,GK/G);
plot_stresses(y,v,we,thetae,Mb,Mt,S_normal,S_shear)

%% TASK5) 
fprintf('<strong> \nTASK5) Divergence and Reversal: \n</strong>')
% Aileron reversal speed + show loss in aileron efficiency with 
% dynamic pressure + decrease of rolling moment coefficient per unit
% aileron deflection with speed (compared to rigid aircraft)

%Do these before, because need to use full Qip
[M,K,Z,Qip,f,CRv,CRd,EI,GK] = otterwing(B,l,b,ts,ba,mhinge,rhop,E,G,nelem,dpm,Ac,hspar);
fprintf('<strong> \nDivergence and Reversal: \n</strong>')

% DIVERGENCE
[udiv,zdiv] = divergence(K,Qip);
fprintf(1,'Divergence speed: %.2f m/s\n', udiv);
if udiv/v_ne > 1.15
    fprintf('Reasonable udiv: (%.2f m/s > 1.15 v_ne)\n',udiv)
end

alpha0 = 1; % deg
alpha = ael_deform(udiv,alpha0,20);
    

% % REVERSAL
[urev,zrev] = reversal(K, Qip, f, CRv, CRd);
fprintf(1,'Reversal speed: %.2f m/s\n', urev);
if urev/v_ne > 1.15
    fprintf('Reasonable urev: (%.2f m/s > 1.15 v_ne)\n',urev)
end

% reversal efficiency analysis 
delta = deg2rad(-1); % given angle of 1 deg of deflection
eff_reversal(urev,S,l,CRd,CRv,delta,K,Qip,f)


%% TASK4) 
fprintf('<strong> \nTASK4) Flutter: \n</strong>')
% Aeroelastic stability analysis, check if design is free of wing flutter
% and static aeroelastic divergenc for velocities up to 1.15*v_ne

% % FLUTTER
% 
% Modal subspace analysis
fprintf('<strong> \nFlutter, Modal subspace analysis: \n</strong>')
% 
[evec, eval] = eig(K,M);
[eval,ind] = sort(diag(eval));
evec = evec(:,ind);
nmode = 10;
Zm = evec(:,1:nmode);
Mz = Zm'*M*Zm;
Kz = Zm'*K*Zm;

ki = length(Qip.ktab);
for i = 1:ki
    Qip.Qtab(1:nmode,1:nmode,i) = Zm'*Qip.Qtab(:,:,i)*Zm;
end
Qip.Qtab = Qip.Qtab(1:nmode,1:nmode,1:ki);

[ucrit, pcrit, zcrit] = flutter(Mz,Kz,Qip,udiv);
fcrit = (ucrit/b)*imag(pcrit)/2/pi;
fprintf(1,'Flutter speed: %.2f m/s \n',ucrit);
fprintf(1,'Frequency of the critical mode: %.2f Hz \n',fcrit);

