function D_el_ip = Tangent_Elastoplastic_DruckerPrager(MESH, indx_p, indx_apex, d_gamma, Etrial, PARAM, nip)

% Compute the consistent tangent modulus for the Drucker Prager Elasto-plastic
% material with piece-wise linear isotropic hardening (plain strain)
% OLD VERSION BELOW derived from HYPLAS2 code (de Souza Neto, E. A.)
% Alban Souche, PGP, UiO, 15/07/2016

% 3D version

nel      = size(MESH.ELEMS,2);
D_el_ip  = zeros(nel,nip,36);

% Initialize entries to elastic tangent modulus
D_el_ip(:,:,1) = bsxfun(@times, 2*PARAM.ep_G(MESH.elem_markers)+PARAM.ep_lambda(MESH.elem_markers), ones(nip,1))' ;
D_el_ip(:,:,2) = bsxfun(@times, PARAM.ep_lambda(MESH.elem_markers), ones(nip,1))';
D_el_ip(:,:,3) = bsxfun(@times, PARAM.ep_lambda(MESH.elem_markers), ones(nip,1))';
D_el_ip(:,:,7) = bsxfun(@times, PARAM.ep_lambda(MESH.elem_markers), ones(nip,1))';
D_el_ip(:,:,8) = bsxfun(@times, 2*PARAM.ep_G(MESH.elem_markers)+PARAM.ep_lambda(MESH.elem_markers), ones(nip,1))' ;
D_el_ip(:,:,9) = bsxfun(@times, PARAM.ep_lambda(MESH.elem_markers), ones(nip,1))';
D_el_ip(:,:,13) = bsxfun(@times, PARAM.ep_lambda(MESH.elem_markers), ones(nip,1))';
D_el_ip(:,:,14) = bsxfun(@times, PARAM.ep_lambda(MESH.elem_markers), ones(nip,1))';
D_el_ip(:,:,15) = bsxfun(@times, 2*PARAM.ep_G(MESH.elem_markers)+PARAM.ep_lambda(MESH.elem_markers), ones(nip,1))' ;
D_el_ip(:,:,22) = bsxfun(@times, PARAM.ep_G(MESH.elem_markers), ones(nip,1))';
D_el_ip(:,:,29) = bsxfun(@times, PARAM.ep_G(MESH.elem_markers), ones(nip,1))';
D_el_ip(:,:,36) = bsxfun(@times, PARAM.ep_G(MESH.elem_markers), ones(nip,1))';

% Index manipulation
indx_cone = indx_p;
indx_cone(indx_apex) = false;
indx_apex = indx_p;
indx_apex(indx_cone) = false;


% Compute the elastoplastic tangent modulus at the itegration point in the plastic domain for the smooth DP cone 
if find(indx_cone==true)
    % Ini. elastoplastic tangent matrix entries
%     soid    = [1 1 0];
%     foid    = [1 0 0;0 1 0; 0 0 .5];
    soid    = [1 1 1 0 0 0];
    foid    = [1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 1 0 0 0 ; 0 0 0 .5 0 0 ; 0 0 0 0 .5 0 ; 0 0 0 0 0 .5 ]; 
    devprj  = foid - (soid'*soid)/3;
    soidmat = soid'*soid;
    R2G     = 2 * PARAM.ep_G;
    ROOT2   = sqrt(2);
    
    
    for ip = 1:nip
        
        % DP cone %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        indxel_ip = indx_cone(ip,:);
        indxel_ip = find(indxel_ip==true); % collecting the relevant elements
        idel_ip   = MESH.elem_markers(indxel_ip);
        
        if indxel_ip
                            
%             Edev      = zeros(4,length(indxel_ip));
%             Eii       = Etrial.xx(ip,indxel_ip)+Etrial.yy(ip,indxel_ip)+Etrial.zz(ip,indxel_ip);
%             Edev(1,:) = Etrial.xx(ip,indxel_ip) - Eii/3;
%             Edev(2,:) = Etrial.yy(ip,indxel_ip) - Eii/3;
%             Edev(3,:) = Etrial.xy(ip,indxel_ip) /2 ; % back to the physical notation
%             Edev(4,:) = Etrial.zz(ip,indxel_ip) - Eii/3;
            
            Edev      = zeros(6,length(indxel_ip));
            Eii       = Etrial.xx(ip,indxel_ip)+Etrial.yy(ip,indxel_ip)+Etrial.zz(ip,indxel_ip);
            Edev(1,:) = Etrial.xx(ip,indxel_ip) - Eii/3;
            Edev(2,:) = Etrial.yy(ip,indxel_ip) - Eii/3;
            Edev(3,:) = Etrial.zz(ip,indxel_ip) - Eii/3;
            Edev(4,:) = Etrial.xy(ip,indxel_ip) /2 ; % back to the physical notation
            Edev(5,:) = Etrial.yz(ip,indxel_ip) /2 ; % back to the physical notation
            Edev(6,:) = Etrial.zx(ip,indxel_ip) /2 ; % back to the physical notation

            d_g    = d_gamma(ip,indxel_ip);
            
%             EdevNORM = sqrt(Edev(1,:).^2 + Edev(2,:).^2 + Edev(4,:).^2 + 2*Edev(3,:).^2);  % ==sqrt(2*(second invariant of strain tensor))
            EdevNORM = sqrt(Edev(1,:).^2 + Edev(2,:).^2 + Edev(3,:).^2 + 2*Edev(4,:).^2 + 2*Edev(5,:).^2 + 2*Edev(6,:).^2);  % ==sqrt(2*(second invariant of strain tensor))
            EdevINV  = 1./EdevNORM;
            EdevINV(EdevNORM==0) = 0;
            
            D = zeros(6,length(indxel_ip));
            D(1,:) = Edev(1,:).*EdevINV;
            D(2,:) = Edev(2,:).*EdevINV;
            D(3,:) = Edev(3,:).*EdevINV;
            D(4,:) = Edev(4,:).*EdevINV;
            D(5,:) = Edev(5,:).*EdevINV;
            D(6,:) = Edev(6,:).*EdevINV;
            
%             Dmat = [ D(1,:).*D(1,:); D(2,:).*D(1,:); D(3,:).*D(1,:); ...
%                      D(1,:).*D(2,:); D(2,:).*D(2,:); D(3,:).*D(2,:); ...
%                      D(1,:).*D(3,:); D(2,:).*D(3,:); D(3,:).*D(3,:)  ]; % Dmat = D*D;        
%             DImat = [ D(1,:)  ; D(2,:)  ; D(3,:)  ; ...
%                       D(1,:)  ; D(2,:)  ; D(3,:)  ; ...
%                       D(1,:)*0; D(2,:)*0; D(3,:)*0]; % DImat = D*I;             
%             IDmat = [ D(1,:); D(1,:); 0*D(1,:); ...
%                       D(2,:); D(2,:); 0*D(2,:); ...
%                       D(3,:); D(3,:); 0*D(3,:) ]; % IDmat = I*D;
                  
            Dmat = [ D(1,:).*D(1,:); D(2,:).*D(1,:); D(3,:).*D(1,:); D(4,:).*D(1,:); D(5,:).*D(1,:); D(6,:).*D(1,:); ...
                     D(1,:).*D(2,:); D(2,:).*D(2,:); D(3,:).*D(2,:); D(4,:).*D(2,:); D(5,:).*D(2,:); D(6,:).*D(2,:); ...
                     D(1,:).*D(3,:); D(2,:).*D(3,:); D(3,:).*D(3,:); D(4,:).*D(3,:); D(5,:).*D(3,:); D(6,:).*D(3,:); ...
                     D(1,:).*D(4,:); D(2,:).*D(4,:); D(3,:).*D(4,:); D(4,:).*D(4,:); D(5,:).*D(4,:); D(6,:).*D(4,:); ...
                     D(1,:).*D(5,:); D(2,:).*D(5,:); D(3,:).*D(5,:); D(4,:).*D(5,:); D(5,:).*D(5,:); D(6,:).*D(5,:); ...
                     D(1,:).*D(6,:); D(2,:).*D(6,:); D(3,:).*D(6,:); D(4,:).*D(6,:); D(5,:).*D(6,:); D(6,:).*D(6,:) ]; % Dmat = D*D;
            
            DImat = [ D(1,:)  ; D(2,:)  ; D(3,:)  ; D(4,:)  ; D(5,:)  ; D(6,:)  ; ...
                      D(1,:)  ; D(2,:)  ; D(3,:)  ; D(4,:)  ; D(5,:)  ; D(6,:)  ; ...
                      D(1,:)  ; D(2,:)  ; D(3,:)  ; D(4,:)  ; D(5,:)  ; D(6,:)  ; ...
                      D(1,:)*0; D(2,:)*0; D(3,:)*0; D(4,:)*0; D(5,:)*0; D(6,:)*0; ...
                      D(1,:)*0; D(2,:)*0; D(3,:)*0; D(4,:)*0; D(5,:)*0; D(6,:)*0; ...
                      D(1,:)*0; D(2,:)*0; D(3,:)*0; D(4,:)*0; D(5,:)*0; D(6,:)*0 ]; % DImat = D*I;      
                         
            IDmat = [ D(1,:); D(1,:); D(1,:); 0*D(1,:); 0*D(1,:); 0*D(1,:); ...
                      D(2,:); D(2,:); D(2,:); 0*D(2,:); 0*D(2,:); 0*D(2,:); ...
                      D(3,:); D(3,:); D(3,:); 0*D(3,:); 0*D(3,:); 0*D(3,:); ...
                      D(4,:); D(4,:); D(4,:); 0*D(4,:); 0*D(4,:); 0*D(4,:); ...
                      D(5,:); D(5,:); D(5,:); 0*D(5,:); 0*D(5,:); 0*D(5,:); ...
                      D(6,:); D(6,:); D(6,:); 0*D(6,:); 0*D(6,:); 0*D(6,:) ]; % IDmat = I*D;
  
            AUX    = 1 ./ (PARAM.ep_G(idel_ip) + PARAM.ep_K(idel_ip)*PARAM.ep_eta*PARAM.ep_etabar + PARAM.ep_H*(PARAM.ep_xi^2));
            AFACT  = R2G(idel_ip) .* ( 1 - ( d_g ./ (ROOT2*EdevNORM) ) );
            BFACT  = R2G(idel_ip) .* ( ( d_g ./ (ROOT2*EdevNORM) ) - (AUX.*PARAM.ep_G(idel_ip)) );
            CFACT  = -ROOT2*PARAM.ep_G(idel_ip).*PARAM.ep_K(idel_ip).*AUX;
            DFACT  = PARAM.ep_K(idel_ip) .* (1 - PARAM.ep_K(idel_ip).*PARAM.ep_eta*PARAM.ep_etabar.*AUX);
            
            tmp1 = devprj(:)*AFACT;
            tmp2 = bsxfun(@times,Dmat,BFACT);
            tmp3 = bsxfun(@times, ( PARAM.ep_eta*DImat + PARAM.ep_etabar*IDmat ), CFACT);
            tmp4 = soidmat(:)*DFACT;
            
            D_el_ip(indxel_ip,ip,:) = (tmp1+tmp2+tmp3+tmp4)';      
            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        % DP apex %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        indxel_ip = indx_apex(ip,:);
        indxel_ip = find(indxel_ip==true); % collecting the relevant elements
        idel_ip   = MESH.elem_markers(indxel_ip);
        
        if indxel_ip
            
            if PARAM.ep_H == 0
                D_el_ip(indxel_ip,ip,:) = 0;
            else
                alpha  = PARAM.ep_xi/PARAM.ep_etabar;
                beta   = PARAM.ep_xi/PARAM.ep_eta;
                AFACT  = PARAM.ep_K(idel_ip).*(1-PARAM.ep_K(idel_ip)./(PARAM.ep_K(idel_ip)+alpha*beta*PARAM.ep_H));
                D_el_ip(indxel_ip,ip,:) = (soidmat(:)*AFACT)';
            end
            
        end       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    end 
end



% %%%%% SLOMO !!!
% 
% nel      = size(MESH.ELEMS,2);
% D_el_ip  = zeros(nel,nip,9);
% 
% % Initialize entries to elastic tangent modulus
% D_el_ip(:,:,1) = 2*PARAM.ep_G+PARAM.ep_lambda;
% D_el_ip(:,:,2) = PARAM.ep_lambda;
% D_el_ip(:,:,4) = PARAM.ep_lambda;
% D_el_ip(:,:,5) = 2*PARAM.ep_G+PARAM.ep_lambda;
% D_el_ip(:,:,9) = PARAM.ep_G;
% 
% indx_tmp = false(size(indx_p));
% indx_tmp(indx_apex) = true;
% all_indx = indx_p + indx_tmp; % 0==elastic 1==DP_cone 2==DP_apex 
% 
% % Ini. elastoplastic tangent matrix entries
% soid   = [1 1 0];
% foid   = [1 0 0;0 1 0; 0 0 .5];
% D_elpl = zeros(3, 3);
% 
% for iel = 1:nel
%     
%     for ip = 1:nip
%       
%         toeva = all_indx(ip,iel);
%         
%         switch toeva
% 
%             case 1 % elastoplastic tangent matrix entries (DP cone)
%                 Ee_ii    = Etrial.xx(ip,iel)+Etrial.yy(ip,iel)+Etrial.zz(ip,iel);
%                 Eedev_xx = Etrial.xx(ip,iel) - Ee_ii/3; 
%                 Eedev_yy = Etrial.yy(ip,iel) - Ee_ii/3; 
%                 Eedev_zz = Etrial.zz(ip,iel) - Ee_ii/3; 
%                 Eedev_xy = Etrial.xy(ip,iel) /2;  % back to physical notation
%                 
%                 Eedev_norm  = sqrt( Eedev_xx^2 + Eedev_yy^2 + Eedev_zz^2 + 2*(Eedev_xy^2) ); % ==sqrt(2*(second invariant of strain tensor))  
%                 Eed_inv = 1/Eedev_norm;
%                 Eed_inv(Eedev_norm==0) = 0; 
%                 
%                 unidev(1) = Eedev_xx*Eed_inv; 
%                 unidev(2) = Eedev_yy*Eed_inv; 
%                 unidev(3) = Eedev_xy*Eed_inv; 
%                 
%                 aux    = 1 / (PARAM.ep_G + PARAM.ep_K*PARAM.ep_eta*PARAM.ep_etabar + PARAM.ep_H*(PARAM.ep_xi^2));
%                 afact  = 2*PARAM.ep_G * (1 - d_gamma(ip,iel)/(sqrt(2)*Eedev_norm));
% %                 afacd3 = afact/3;
%                 bfact  = 2*PARAM.ep_G * (d_gamma(ip,iel)/(sqrt(2)*Eedev_norm) - PARAM.ep_G*aux);
%                 cfact  = -sqrt(2) * PARAM.ep_G * PARAM.ep_K * aux;
%                 dfact  = PARAM.ep_K * (1 - PARAM.ep_K*PARAM.ep_eta*PARAM.ep_etabar*aux);
%                 
%                 for i = 1:3
%                     for j=1:3
%                         % D_elpl(i,j) = afact*foid(i,j) + bfact*unidev(i)*unidev(j) ...
%                         %             + cfact*( PARAM.ep_nu_phi*unidev(i)*soid(j) +  PARAM.ep_nu_dil*soid(i)*unidev(j)) ...
%                         %             + (dfact-afacd3)*soid(i)*soid(j);
%                         
%                         devprj = foid(i,j) - ( soid(i)*soid(j) /3 ); 
%                         D_elpl(i,j) = afact*devprj  + bfact*unidev(i)*unidev(j) ...
%                             + cfact*( PARAM.ep_eta*unidev(i)*soid(j) +  PARAM.ep_etabar*soid(i)*unidev(j) ) ...
%                             + dfact*soid(i)*soid(j);
%                     end
%                 end
%                 
%                 D_el_ip(iel,ip,:) = D_elpl(:);  
%                 
%                 
%             case 2 % elastoplastic tangent matrix entries (DP apex)
%                 alpha  = PARAM.ep_xi/PARAM.ep_etabar;
%                 beta   = PARAM.ep_xi/PARAM.ep_eta;
%                 afact  = PARAM.ep_K*(1-PARAM.ep_K/(PARAM.ep_K+alpha*beta*PARAM.ep_H));
% %                 D_elpl = afact * ((soid')*soid);    
%                 for i = 1:3
%                     for j=1:3
%                         D_elpl(i,j) = afact*soid(i)*soid(j);
%                     end
%                 end
%                 
%                 D_el_ip(iel,ip,:) = D_elpl(:);   
%         end
%         
%    
%     end
% end

