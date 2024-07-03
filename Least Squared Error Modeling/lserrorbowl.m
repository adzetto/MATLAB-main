% 2-coefficient error surface
Ho=[-1.5 -2.5]';	% actual filter

[HX, HY]=meshgrid(-6.5:.2:+3.5, -7.5:.2:2.5);

err = 15.*ones(size(HX)) - 2.*(HX.*Ho(1) + HY.*Ho(2)) + HX.^2 + HY.^2;
errmin = 15 - Ho(1)^2 - Ho(2)^2;

figure(1);

meshc(HX,HY, err);

ax3=axis;
zmin=ax3(5);

grid, xlabel('h0'), ylabel('h1'), zlabel('squared error');
title('Squared-Error Surface for 2-Parameter System');
hold on;
plot3([Ho(1) Ho(1)],[Ho(2) Ho(2)],[zmin errmin],'k:');
hold off;
view(-30,10);
caxis([0 1]);
grid on;
