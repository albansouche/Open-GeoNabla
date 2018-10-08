function d_gamma = Return_mapping_VonMises(MESH, F, param,  indx_p)

% Part of Open-GeoNabla, copyright GPLv3, 2018
% https://github.com/albansouche/Open-GeoNabla/
% Physics of Geological Processes (PGP) , The NJORD Centre, Dept of Geosciences, University of Oslo
% Author: Alban Souche

% Return mapping algorithm for Von Mises material with piece-wise linear
% isotropic hardening (plain strain)
% Inspired from HYPLAS2 code (de Souza Neto, E. A.) 


indx_el_p = repmat(MESH.elem_markers,size(indx_p,1),1);
indx_el_p = indx_el_p(indx_p);
indx_p    = find(indx_p);

% Return mapping Von Mises
numer_tmp = F(indx_p);
G         = param.ep_G(indx_el_p);  G  = G(:);
denom_tmp = -3*G - param.ep_H;
d_gamma   = - numer_tmp./denom_tmp ;




