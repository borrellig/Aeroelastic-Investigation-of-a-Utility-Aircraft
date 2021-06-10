function  eff_reversal(urev,S,l,CRd,CRv,delta,K,Qip,f)
% possible to check efficiency with other formula

rho = Qip.rho;
u = linspace(0,urev,20);
q = 0.5*rho*u.^2;
Q0 = Qip.Qtab(:,:,1);
f0 = f;

for i = 1:length(q)
    % displacement for a given velocity (delta is given)
    v(:,i) = (K-q(i)*Q0)\(q(i)*f0*delta);
    
    % Rolling moment due to rigid body motion
    R_rigid(i) = q(i)*S*l*CRd*delta;
    
    % Rolling moment accounting also for flexible structure
    R_flex(i) = q(i)*S*l*CRv*v(:,i)+R_rigid(i);
    Cr_flex(i) = R_flex(i)/(q(i)*S*l);
end

eta = R_flex./R_rigid;

% plot efficiency
figure
%subplot(2,1,1)
plot(q,eta,'o-')
title('Rolling efficiency','FontSize',16)
xlabel('q_{dyn} [N/m^2]','FontSize',15)
ylabel('Efficiency','FontSize',15)

% subplot(3,1,2)
% title('Rolling coefficient')
% xlabel('u')
% ylabel('C_{R}')
% hold on
% plot(u,Cr_flex,'o-')
% grid on

%subplot(2,1,2)
figure
plot(u,R_flex,'o-',u,R_rigid,'--')
title('Rolling Moment','FontSize',16)
xlabel('u [m/s]','FontSize',15)
ylabel('Rolling Moment [Nm]','FontSize',15)
legend('Flexible wing','Rigid wing','FontSize',13)

