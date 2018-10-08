function alpha = poly_search(f,res0)
%function alpha = golden_section_search2(f,beta)

% Find residuum for s=0, s=1/2, s=1
% f0      = f(0); 
% f1      = f(1);
% fmid    = f(1/2);
% Find coefficients in the equation: a*x^2+b*x+c=0
% c   = f0;
% tmp = [1 1; 1/4 1/2]\[f1-f0; fmid-f0];
% a   = tmp(1);
% b   = tmp(2);

% find residum for s=0, s=1/3, s=2/3, s=1
% f03     = f(0/3); 
% f13     = f(1/3);
% f23     = f(2/3);
% f33     = f(3/3);

% f00 = f(0);
f00 = res0;
% f01 = f(01/10); 
f02 = f(02/10);
% f03 = f(03/10);
f04 = f(04/10);
% f05 = f(05/10);
f06 = f(06/10);
% f07 = f(07/10);
f08 = f(08/10);
% f09 = f(09/10);
f10 = f(10/10);

% Find coefficients in the equation using four points: a*x^4+b*x^3+c*x^2+d*x+e=0
e2  = f00;
tmp = [(0.2)^4 (0.2)^3 (0.2)^2 0.2;...
       (0.4)^4 (0.4)^3 (0.4)^2 0.4;...
       (0.6)^4 (0.6)^3 (0.6)^2 0.6;...
       (0.8)^4 (0.8)^3 (0.8)^2 0.8;...
           1^4     1^3     1^2   1]\ ...
      [f02-f00; f04-f00; f06-f00; f08-f00; f10-f00];
a2  = tmp(1);
b2  = tmp(2);
c2  = tmp(3);
d2  = tmp(4);

xx = linspace(0,1,1e4);
[~, idx] = min(a2*xx.^4+b2*xx.^3+c2*xx.^2+d2*xx+e2);
x_ext4   = xx(idx);

% figure(22)
% clf; hold on
% % plot(0:0.1:1,[f00 f01 f02 f03 f04 f05 f06 f07 f08 f09 f10],'b.');
% % plot(0:0.5:1,[f0 fmid f1],'ob')
% plot(0:0.2:1,[f00 f02 f04 f06 f08 f10],'go');
% % plot(xx,a*xx.^2+b*xx+c,'r-')
% plot(xx,a2*xx.^4+b2*xx.^3+c2*xx.^2+d2*xx+e2,'g-')
% % Find exteremum
% % x_ext  = -b/(2*a);
% % plot(x_ext, f(x_ext), '*r')
% plot(x_ext4, f(x_ext4), '*r')
% drawnow
% % plot(beta,f(beta),'ks')
% xlabel('alpha')
% ylabel('f')
% box on
% alpha_new = x_ext4
if abs(x_ext4-0.5)<0.5 % minimum in the range of (0,1)
    alpha = x_ext4;
else % minimum at the edge
    if f00<f10 
        alpha = 0;
    else 
        alpha = 1;
    end
end