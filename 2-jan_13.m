% 3-D Spiral
clc
clear

x = 0:0.1:10*pi;
y1 = sin(x);
y2 = cos(x);
y3=x;
figure
plot3(y1, y2,y3,'r--','LineWidth',2)

xlabel('X-axis','FontSize',10,'Color','#00FFFF')
ylabel('Y-axis','FontSize',20,'Color','#0F00F0')
zlabel('Z-axis')

title('3D Curve: sin(x) vs cos(x)','FontSize',25,'Color','r','BackgroundColor','#000000')

grid on
legend('x, sin(x), cos(x)','Location','best')
view(45,30)
