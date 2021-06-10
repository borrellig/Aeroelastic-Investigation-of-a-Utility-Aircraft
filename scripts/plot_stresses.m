function plot_stresses(y,v,we,thetae,Mb,Mt,S_normal,S_shear)

%% plots to compare real with approx displacements
% approximation is always larger wrt real one (good because conservative-?)
% 
% considering that stresses is defined inside the element, we cannot plot
% the last element with finite differences. This explains the ranges taken
% into account for the plots (basically to have symmetry)

w = v(1:3:end-2);
theta = v(3:3:end);

figure
subplot(2,1,1)
title('w [m]')
hold on
plot(y(2:end-1),we(2:end))
plot(y(2:end),w(2:end),'--')
grid on
legend('with shape function','no shape function')

subplot(2,1,2)
title('theta [rad]')
hold on
plot(y(2:end-1),thetae(2:end))
plot(y(2:end),theta(2:end),'--')
grid on

%% Plot bending and twist moments
figure
subplot(2,1,1)
title('Mb [Nm]')
hold on
plot(y(2:end-1),Mb(2:end))
grid on

subplot(2,1,2)
title('Mt [Nm]')
hold on
plot(y(2:end-1),Mt(2:end))
grid on

%% Plot Normal and Shear stresses for different locations
figure
subplot(2,1,1)
title('Normal stresses [N/m^2]')
hold on
% plot(y(1:end-1),S_normal(1,:),'*-')
plot(y(2:end-1),S_normal(:,2:end))
grid on

subplot(2,1,2)
title('Shear stresses [N/m^2]')
hold on
% plot(y(1:end-1),S_shear(1,:),'*-')
plot(y(2:end-1),S_shear(:,2:end))
grid on

% Plot most critical stresses (tension) [if plotting compression -> buckling]
% in tension stresses are positive
figure
subplot(2,1,1)
plot(y(2:end-1),S_normal(1,2:end),'k','LineWidth',2)
title('Normal stresses distribution','fontsize',15)
ylabel('$\sigma_N$ [N/m$^2$]','Interpreter','latex','Fontsize',15)
xlabel('y [m]','fontsize',15)
% plot(y(1:end-1),S_normal)


subplot(2,1,2)
plot(y(2:end-1),S_shear(1,2:end),'r','LineWidth',2)
title('Shear stresses distribution','fontsize',15)
ylabel('$\sigma_S$ [N/m$^2$]','Interpreter','latex','Fontsize',15)
xlabel('y [m]','fontsize',15)
% plot(y(1:end-1),S_shear)
