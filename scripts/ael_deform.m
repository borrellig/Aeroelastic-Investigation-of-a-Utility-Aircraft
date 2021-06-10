function alpha = ael_deform(udiv,alpha0,n)
% function to compute the deformation, namely the alpha_elastic, which is
% the quantity of interest in aeroelasticity

rho = 1.225;
qdiv = 0.5*rho*udiv^2;
q = linspace(0,qdiv,n);
u = linspace(0,udiv,n);

for i = 1:n
    alpha(i) = 1/(qdiv/q(i)-1)*alpha0;
end

figure
plot(u,alpha,'ko-',u,alpha,'b',udiv*ones(n),alpha(1:n),'r--')
xlabel('u [m/s]', 'fontsize',15)
ylabel('$\alpha$ [deg]','Interpreter','latex','Fontsize',15)

% when q approaches qdiv the deformation goes to infinite as expected

