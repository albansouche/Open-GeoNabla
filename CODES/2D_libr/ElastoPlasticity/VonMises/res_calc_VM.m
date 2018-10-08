function res = res_calc_VM(MESH, InvJ_all, DetJ_all, Etot, param, ieuv, nip, nelblo, BC, indx_p, dU, Fext)

% Part of Open-GeoNabla, copyright GPLv3, 2018
% https://github.com/albansouche/Open-GeoNabla/
% Physics of Geological Processes (PGP) , The NJORD Centre, Dept of Geosciences, University of Oslo
% Author: Alban Souche


nel = size(MESH.ELEMS,2);

dEtrial      = strain_tot_calc(MESH, InvJ_all, dU, ieuv, nelblo);
dEtrial.xy   = 2*dEtrial.xy; % gamma convention
% Add strain increment
Etrial.xx = Etot.xx + dEtrial.xx;
Etrial.yy = Etot.yy + dEtrial.yy;
Etrial.zz = Etot.zz + dEtrial.zz;
Etrial.xy = Etot.xy + dEtrial.xy;
Etrial.ii = Etrial.xx+Etrial.yy+Etrial.zz;
% Calculate Pressure
P       = param.K*Etrial.ii;
% Calculate deviatoric stress
Sdev.xx = 2*param.G * ( Etrial.xx - Etrial.ii/3 );
Sdev.yy = 2*param.G * ( Etrial.yy - Etrial.ii/3 );
Sdev.zz = 2*param.G * ( Etrial.zz - Etrial.ii/3 );
Sdev.xy = 2*param.G * ( Etrial.xy/2 );

% Calculate J2
J2 = (Sdev.xx.^2 + Sdev.yy.^2 + Sdev.zz.^2)/2 + Sdev.xy.^2;
% Calculate Von Mises effective stress
Q  = sqrt(3*J2);

% Check for plasticity consistency %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F      = Q - param.S_y;
% indx_p = F>0;

d_gamma = zeros(nip,nel);

if find(indx_p)
    % Return mapping
    d_g = Return_mapping_VonMises(F, param, indx_p);
    d_gamma(indx_p) = d_g;
end

% Apply correction to the stress tensor
fact = 1 - 3*param.G*d_gamma./Q;
Sdev.xx = fact.*Sdev.xx;
Sdev.yy = fact.*Sdev.yy;
Sdev.zz = fact.*Sdev.zz;
Sdev.xy = fact.*Sdev.xy;

Fint = internal_force(MESH, InvJ_all, DetJ_all, Sdev, P, nip, nelblo);

res = Fext-Fint;
res(BC.ind)  = 0;
