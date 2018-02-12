function [ipx, ipw] = ip_tetra(nip)

%%% BASED ON IP_TRIANGLE.m from MILAMIN
% Keast, P., 1986, Moderate Degree Tetrahedral Quadrature Formulas, 
% Comm. Appl. Num. Meth., 55, pp. 339???348

switch nip
    case 1
        ipx(1,1) = 1/4;
        ipx(1,2) = 1/4;
        ipx(1,3) = 1/4;

        ipw(1)   = 1/6;

    case 4
        a = 0.5854101966249685; %(5 + 3*sqrt(5))/20;
        b = 0.1381966011250105; %(5 - sqrt(5))/20; 
 
        ipx(1,1:3) = [a ,b ,b];
        ipx(2,1:3) = [b ,a b];
        ipx(3,1:3) = [b ,b ,a];
        ipx(4,1:3) = [b ,b ,b];

        ipw(1)  = 1/24;
        ipw(2)  = 1/24;
        ipw(3)  = 1/24;
        ipw(4)  = 1/24;

    case 5
        ipx(1,1:3) = [1/4 ,1/4 ,1/4];
        ipx(2,1:3) = [1/2 ,1/6 ,1/6];
        ipx(3,1:3) = [1/6 ,1/2 ,1/6];
        ipx(4,1:3) = [1/6 ,1/6 ,1/2];
        ipx(5,1:3) = [1/6 ,1/6 ,1/6];
        
        ipw(1)  = -4/30;
        ipw(2)  = 9/120;
        ipw(3)  = 9/120;
        ipw(4)  = 9/120;
        ipw(5)  = 9/120;
        
    case 11
        a = 0.3994035761667992; % (1 + sqrt(5/14))/4; 
        b = 0.1005964238332008; % (1 - sqrt(5/14))/4; 
        ipx(1,1:3) = [1/4   ,1/4   ,1/4  ];
        ipx(2,1:3) = [11/14 ,1/14  ,1/14 ];
        ipx(3,1:3) = [1/14  ,11/14 ,1/14 ];
        ipx(4,1:3) = [1/14  ,1/14  ,11/14];
        ipx(5,1:3) = [1/14  ,1/14  ,1/14 ];
        ipx(6,1:3) = [a ,a ,b];
        ipx(7,1:3) = [a ,b ,a];
        ipx(8,1:3) = [a ,b ,b];
        ipx(9,1:3) = [b ,a ,a];
        ipx(10,1:3) = [b ,a ,b];
        ipx(11,1:3) = [b ,b ,a];
        
        ipw(1)  = -74/5625;
        ipw(2)  = 343/45000;
        ipw(3)  = 343/45000;
        ipw(4)  = 343/45000;
        ipw(5)  = 343/45000;
        ipw(6)  = 56/2250;
        ipw(7)  = 56/2250;
        ipw(8)  = 56/2250;
        ipw(9)  = 56/2250;
        ipw(10) = 56/2250;
        ipw(11) = 56/2250;
      
    otherwise
        error('Unknown integration rule')
        
end

end