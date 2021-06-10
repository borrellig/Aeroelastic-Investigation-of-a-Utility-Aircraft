function [ucrit,pcrit,zcrit] = flutter(M,K,Qip,udiv)
% [ucrit,pcrit,zcrit] = flutter(M,K,Qip)
%
% Compute flutter speed, frequency and critical mode.
%
% M : mass matrix
% K : stiffness matrix
% Qip : aerodynamic loads interpolation struct
%
% (c) 2004-2016 Dan Borglund <dodde@kth.se> and David Eller <dlr@kth.se>

  
  % set return values
  
  umin = 0.1;
  u = linspace(umin,udiv,40);
  
for i = 1:length(u)
      
      [kbounds] = pk_bounds(u(i),M,K,Qip,10);
      for j = 1:length(kbounds)-1
          [pu,vu] = pk_bisect(j,u(i),M,K,Qip,kbounds);
          pp(j) = pu;
          vv(:,j) = vu;
      end
      
      p(i,:) = pp;
      v(i,:,:) = vv;
     
end

ind = min(find(real(p) > 0))-1;
L = length(u);
col = ceil(ind/L);
row = mod(ind,L);
pcrit = p(row,col);
% fprintf('Critical eigenvalue related to \nreduced freq. %d\n',col)

ucrit = u(row);
zcrit = v(row,:,col);

% figure
% xlabel('Re(p)')
% ylabel('Im(p)')
% hold on
% plot(real(p(:,1)),imag(p(:,1)),'b*',real(p(1,1)),imag(p(1,1)),'ro',real(p(end,end)),imag(p(end,end)),'ro')
% plot(real(p(:,2)),imag(p(:,2)),'r*')
% plot(real(p(:,3)),imag(p(:,3)),'g*')
% plot(real(p(:,4)),imag(p(:,4)),'k*')
% plot(real(p(:,5)),imag(p(:,5)),'m*')
% legend('k1','k2','k3','k4','k5')
% grid on

figure
plot(u,real(p(:,1)),'bo',u(1:end-5),real(p(1:end-5,2)),'ro',u,real(p(:,3)),'o',u,real(p(:,4)),'ko')
% plot(u,real(p(:,5)),'m*')
% plot(u,real(p(:,6)),'c*')
% plot([u(1),u(end)],[0 0],'k')
legend('k1','k2','k3','k4','FontSize',14)
xlabel('Speed [m/s]','FontSize',15)
ylabel('Re(p)','FontSize',15)