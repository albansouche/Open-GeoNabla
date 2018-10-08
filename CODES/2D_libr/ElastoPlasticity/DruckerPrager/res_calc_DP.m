function res = res_calc_DP(MESH, InvJ_all, DetJ_all, Etot, Epl, PARAM, ieuv, BC, indx_p, dU, Fext, P_litho, SOLVER)

% Part of Open-GeoNabla, copyright GPLv3, 2018
% https://github.com/albansouche/Open-GeoNabla/
% Physics of Geological Processes (PGP) , The NJORD Centre, Dept of Geosciences, University of Oslo
% Author: Alban Souche


dEtrial      = strain_tot_calc(MESH, InvJ_all, dU, ieuv, SOLVER.nelblo);
dEtrial.xy   = 2*dEtrial.xy; % gamma convention
% Add strain increment
Etrial.xx = Etot.xx + dEtrial.xx;
Etrial.yy = Etot.yy + dEtrial.yy;
Etrial.zz = Etot.zz + dEtrial.zz;
Etrial.xy = Etot.xy + dEtrial.xy;
Etrial.ii = Etrial.xx+Etrial.yy+Etrial.zz;
% Calculate Pressure
P = bsxfun(@times,PARAM.ep_K(MESH.elem_markers), Etrial.ii) - P_litho;  %% add (-P_litho) negative because of convention
% Calculate deviatoric stress
Sdev.xx = bsxfun(@times , 2*PARAM.ep_G(MESH.elem_markers) , ( Etrial.xx - Etrial.ii/3 )) ;
Sdev.yy = bsxfun(@times , 2*PARAM.ep_G(MESH.elem_markers) , ( Etrial.yy - Etrial.ii/3 )) ;
Sdev.zz = bsxfun(@times , 2*PARAM.ep_G(MESH.elem_markers) , ( Etrial.zz - Etrial.ii/3 )) ;
Sdev.xy = bsxfun(@times , 2*PARAM.ep_G(MESH.elem_markers) , ( Etrial.xy/2 )) ;

% Calculate J2
J2 = (Sdev.xx.^2 + Sdev.yy.^2 + Sdev.zz.^2)/2 + Sdev.xy.^2;

% Check for plasticity consistency %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cohe   = bsxfun(@plus, PARAM.ep_Co(MESH.elem_markers), PARAM.ep_H*Epl.bar_trial);
F      = sqrt(J2) - PARAM.ep_eta*P - PARAM.ep_xi.*cohe;

% Return mapping with update of deviatoric stress tensor and P
% [~, Sdev, P, ~, ~] = Return_mapping_DruckerPrager(MESH, F, PARAM, indx_p, J2, Sdev, P, Epl);
[~, Sdev, P, ~, ~] = Return_mapping_DruckerPrager_cohesion_noise(MESH, F, PARAM, indx_p, J2, Sdev, P, Epl);


Fint = internal_force(MESH, InvJ_all, DetJ_all, Sdev, P + P_litho, SOLVER);

res = Fext - Fint;
res(BC.ind)  = 0;
