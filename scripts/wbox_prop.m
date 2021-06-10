function [l_wbox,xle_wbx,xte_wbx,le_wbx,te_wbx,x_ud,l_up,l_down] = wbox_prop(b)

% function to calculate the properites of the wingbox
% NACA 63A-516, wing profile
load naca63a516.dat
load wingbox.txt
% figure
% axis('equal')
% hold on
% plot(naca63a516(:,1),naca63a516(:,2))
% plot(wingbox(:,1),wingbox(:,2))
% set(gca,'XColor', 'none','YColor','none')

% obtain dimensions of wingbox (contour), to evaluate the area of the wingbox
x = wingbox(:,1);
y = wingbox(:,2);

for i = 1:length(x)-1
    d(i) = sqrt((x(i+1)-x(i))^2+(y(i+1)-y(i))^2);
end
l_wbox = sum(d)*2*b;  % = 2.9458
xle_wbx = (0.1-0.5)*2*b;
xte_wbx = (0.75-0.5)*2*b;
le_wbx = 0.1065*2*b;
te_wbx = 0.0795*2*b;
l_up = sum(d(1:16))*2*b;
l_down = sum(d(18:end-1))*2*b;
x_ud = 0.5*(xle_wbx+xte_wbx);
% check = l_up+l_down+le_wbx+te_wbx %check of length of wingbox

