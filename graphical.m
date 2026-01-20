clc;
clear;
A = [-1 -3; 1 1; 1 -1];
B = [10; 6; 2];
y1 = 0 : 1 : max(B);
Y = y1;

X11 = (B(1)-A(1,1).*Y)./A(1,2);
X12 = (B(2)-A(2,1).*Y)./A(2,2);
X13 = (B(3)-A(3,1).*Y)./A(3,2);

X11 = max(0, X11);
X12 = max(0, X12);
X13 = max(0, X13);

plot(y1, X11, 'b-', y1, X12, 'r--', y1, X13, 'g:')

% corner point
% cx1 = find(Y==0);
% c1 = find(X11 == 0);
% disp(cx1)
% disp(c1)
% line1 = [Y(:,[c1,cx1]);x11(:,[c1,cx1])]';

% c2 = find(X11 == 0);
% line2 = [Y(:,[c2,cx1]);x21(:,[c2,cx1])]';

% c3 = find(X11 == 0);
% line3 = [Y(:,[c3,cx1]);x31(:,[c3,cx1])]';




% corner point
cx1 = find(Y==0);
c1 = find(X11 == 0);
c2 = find(X12 == 0);
c3 = find(X13 == 0);

 line1 = [Y([c1,cx1]); X11([c1,cx1])]';
 line2 = [Y([c2,cx1]); X12([c2,cx1])]';
 line3 = [Y([c3,cx1]); X13([c3,cx1])]';

corpt = unique([line1; line2; line3],'rows');
disp(corpt)
