function [d_gamma, Sdev, P, Epl, indx_apex] = Return_mapping_DruckerPrager_cohesion_noise(MESH, F, PARAM, indx_p, J2, Sdev, P, Epl)

% Part of Open-GeoNabla, copyright GPLv3, 2018
% https://github.com/albansouche/Open-GeoNabla/
% Physics of Geological Processes (PGP) , The NJORD Centre, Dept of Geosciences, University of Oslo
% Author: Alban Souche

% Return mapping algorithm for Drucker Prager material with piece-wise linear
% isotropic hardening (plain strain)
% Inspired from HYPLAS2 code (de Souza Neto, E. A.) 


co_noise = repmat(PARAM.noise,size(indx_p,1),1);
% co_noise = co_noise(indx_p)';
co_noise = co_noise(indx_p);
co_noise = co_noise(:)';

indx_el_p = repmat(MESH.elem_markers,size(indx_p,1),1);
indx_el_p = indx_el_p(indx_p);
indx_p    = find(indx_p);

indx_el_p = indx_el_p(:)';
indx_p = indx_p(:)';

% Return mapping Drucker-Prager on smooth cone
numer_tmp = F(indx_p);  numer_tmp = numer_tmp(:)';
G         = PARAM.ep_G(indx_el_p);  G  = G(:)';
K         = PARAM.ep_K(indx_el_p);  K  = K(:)';
Co        = PARAM.ep_Co(indx_el_p) + co_noise ; Co = Co(:)';
denom_tmp = PARAM.ep_H*PARAM.ep_xi^2 + G + K.*PARAM.ep_eta.*PARAM.ep_etabar;
d_gamma   = numer_tmp./denom_tmp ;

if PARAM.ep_xsi==0
%     No return mapping to the apex needed
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
fact(J2(indx_cone)==0) = 1; % make sure that there is no "inf" value
fact(indx_apex_tmp) = 0; % set zero at the apex

% Update deviatoric stress tensor
Sdev.xx(indx_p)  = fact.*Sdev.xx(indx_p);
Sdev.yy(indx_p)  = fact.*Sdev.yy(indx_p);
Sdev.zz(indx_p)  = fact.*Sdev.zz(indx_p);
Sdev.xy(indx_p)  = fact.*Sdev.xy(indx_p);

% figure(1111),clf, hold on
% 
% for ii= 1:100
% 
% pts_indx = indx_p(ii);
% 
% Sxx = Sdev.xx(pts_indx) + P(pts_indx);
% Syy = Sdev.yy(pts_indx) + P(pts_indx);
% Sxy = Sdev.xy(pts_indx);
%  
% S1    = (Sxx+Syy)/2 - sqrt( ( (Sxx-Syy)^2)/4 + Sxy^2 );
% S2    = (Sxx+Syy)/2 + sqrt( ( (Sxx-Syy)^2)/4 + Sxy^2 );
% 
% Savg  = (S1+S2)/2;
% R     = abs(S1 - Savg);
% 
% % Plot stresses
% plot(S1, 0, 'rd')
% plot(S2, 0, 'bd')
% %     plot(P(pts_indx), 0, 'sk' )
% plot(Savg, 0, 'sc' )
% %     plot( Savg+(R*cosd(2*fric_angle)) , Tau, 'gd')
% %Plot circle
% theta = linspace(0, pi, 100);
% plot( R*cos(theta)+Savg , R*sin(theta) ,'-k' )
% axis equal, grid on
% set(gca,'xdir','reverse')
%     
% x_line = linspace(-3e7,Co(ii)/tand(30),100);
% plot(x_line, -x_line*tand(30)+ Co(ii),'-r')
%     
% 
% drawnow
% end




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


