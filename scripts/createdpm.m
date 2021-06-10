function dpm = createdpm(b)

% creation of matrix with concentrated masses
% reference with respect to the half chord, at this stage

npmass = 5;
dpm = zeros(npmass,5); 

% engine nacelles (gas turbines + propellers)
dJ_eng = [250, 0; 0, 650];

dpm(1,:) = [330, -b-1.3, 2.85, dJ_eng(1,1), dJ_eng(2,2)];
dpm(2,:) = [330, -b-1.3, -2.85, dJ_eng(1,1), dJ_eng(2,2)]; % !!! comment if
% test with clamped

% tail + fuselage (try to separate contribution from tail and fuselage)
m_ft = 1710;
x_ft = -b-0.1;
r_gx = 0.5; % radius of gyration
r_gy = 1.7; 
% dJ_ft = m_ft*[r_gx^2, 0; 0, r_gy^2];
% dpm(3,:) = [m_ft, -b-0.1, 0, dJ_ft(1,1), dJ_ft(2,2)];
% tail (same material of main wing)
rhop = 2780;   S_t = 9.3;    
t_t = 0.01; % assuming thickness = 10mm (spar + skin)
m_tail = rhop*S_t*t_t;
x_tail = -b+8.5;
dpm(3,:) = [m_tail, x_tail, 0, m_tail*r_gx^2, m_tail*r_gy^2];
% fuselage
m_fus = m_ft-m_tail;
x_fus = (m_ft*x_ft-m_tail*x_tail)/m_fus;
dpm(4,:) = [m_fus, x_fus, 0, m_fus*r_gx^2, m_fus*r_gy^2];


% possible to accound for same radii of gyrations of tail+fuselage for the
% payload and the fuel with center of mass in front of wing LE of 10 cm

% payload (minimum of 100kg)
pload = 100;
dJ_pload = pload*[r_gx^2, 0; 0, r_gy^2];
dpm(5,:) = [pload, -b-0.1, 0, dJ_pload(1,1), dJ_pload(2,2)];

% fuel
mfuel = 1170;
dJ_fuel = mfuel*[r_gx^2, 0; 0, r_gy^2];
dpm(6,:) = [mfuel, -b-0.1, 0 , dJ_fuel(1,1), dJ_fuel(2,2)];