function [d_gamma, Sdev, P, Epl, indx_apex] = Return_mapping_DruckerPrager(MESH, F, PARAM, indx_p, J2, Sdev, P, Epl)

% Part of OpenGeoNabla, copyright GPLv3, 2018
% https://github.com/albansouche/OpenGeoNabla/
% Physics of Geological Processes (PGP) , The NJORD Centre, Dept of Geosciences, University of Oslo
% Author: Alban Souche

% Return mapping algorithm for Drucker Prager material with piece-wise linear
% isotropic hardening (plain strain)
% Inspired from HYPLAS2_v2.0 code (de Souza Neto, E. A.) 

indx_el_p = repmat(MESH.elem_markers,size(indx_p,1),1);
indx_p    = find(indx_p);
indx_el_p = indx_el_p(indx_p);

numer_tmp = F(indx_p);

G         = PARAM.ep_G(indx_el_p);     
K         = PARAM.ep_K(indx_el_p);    
Co        = PARAM.ep_Co(indx_el_p);   
denom_tmp = PARAM.ep_H*PARAM.ep_xi^2 + G + K.*PARAM.ep_eta.*PARAM.ep_etabar;

% fix size array for tet4 and tet10 element (comes from the size of F(indx_p))
G = reshape(G , size(numer_tmp)); 
K = reshape(K , size(numer_tmp));
Co = reshape(Co , size(numer_tmp));
denom_tmp = reshape(denom_tmp , size(numer_tmp));

d_gamma   = numer_tmp./denom_tmp ;


if PARAM.ep_xsi==0
    % No return mapping to the apex needed
    indx_apex_tmp = false(size(indx_p));
    indx_apex     = indx_p(indx_apex_tmp);
    indx_cone     = indx_p(~indx_apex_tmp);
else
    % find index on smooth cone and apex
    % Check invalidity of the return mapping on the smooth cone
    indx_apex_tmp = (sqrt(J2(indx_p)) - G.*d_gamma) < 0 ;
    indx_apex     = indx_p(indx_apex_tmp);
    indx_cone     = indx_p(~indx_apex_tmp);
end


% Update P on the smooth cone
P(indx_cone)      = P(indx_cone) - K(~indx_apex_tmp).*PARAM.ep_etabar.*d_gamma(~indx_apex_tmp);
% Update integrated plastic strain on the smooth cone
Epl.bar_trial(indx_cone) = Epl.bar(indx_cone) + PARAM.ep_xi*d_gamma(~indx_apex_tmp);

% Factor for stresses update (cone & apex with zeros at apex)
fact = 1 - ( G .* d_gamma ./ sqrt(J2(indx_p)) );
fact(J2(indx_cone)==0) = 0; % make sure that there is no "inf" value
fact(indx_apex_tmp) = 0; % set zero at the apex

% Update deviatoric stress tensor
Sdev.xx(indx_p)  = fact.*Sdev.xx(indx_p);
Sdev.yy(indx_p)  = fact.*Sdev.yy(indx_p);
Sdev.zz(indx_p)  = fact.*Sdev.zz(indx_p);
Sdev.xy(indx_p)  = fact.*Sdev.xy(indx_p);
Sdev.yz(indx_p)  = fact.*Sdev.yz(indx_p);
Sdev.zx(indx_p)  = fact.*Sdev.zx(indx_p);


% Return mapping Drucker-Prager on apex
if ~isempty(indx_apex)
    
    alpha     = PARAM.ep_xi/PARAM.ep_etabar;
    beta      = PARAM.ep_xi/PARAM.ep_eta;
 
    numer_tmp = - beta*(Co(indx_apex_tmp) + PARAM.ep_H*Epl.bar(indx_apex)) + P(indx_apex);
    denom_tmp = beta.*PARAM.ep_H.*alpha + K(indx_apex_tmp);
    d_Epv     = numer_tmp./denom_tmp;
    
    % Updates at the apex
    P(indx_apex)             = P(indx_apex) - K(indx_apex_tmp).*d_Epv;
    
    Epl.bar_trial(indx_apex) = Epl.bar(indx_apex) + alpha*d_Epv;
    d_gamma(indx_apex_tmp)   = d_Epv/PARAM.ep_etabar; 
    
    check = d_gamma(indx_apex_tmp);
    if check<0
        display('Return Mapping Failed: probably because of too large strain increment');
    end
    
end


