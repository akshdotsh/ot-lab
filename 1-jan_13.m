clc
clear

x = 0:0.1:10;
y1 = sin(x);
y2 = cos(x);

plot(x, y1,'r--')
hold on
plot(x, y2,'g--')
hold off

xlabel('x-axis','FontSize',10,'Color','#00FFFF')
ylabel('y-axis','FontSize',20,'Color','#0F00F0')
title('sinx vs cosx','FontSize',25,'Color','r','BackgroundColor','#000000')

grid on
legend('y = sin(x)','y = cos(x)','Location','best')
