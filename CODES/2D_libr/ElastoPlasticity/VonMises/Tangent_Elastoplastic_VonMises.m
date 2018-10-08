function D_el_ip = Tangent_Elastoplastic_VonMises(MESH, indx_p, d_gamma, Sdev, param, nip)

% Part of Open-GeoNabla, copyright GPLv3, 2018
% https://github.com/albansouche/Open-GeoNabla/
% Physics of Geological Processes (PGP) , The NJORD Centre, Dept of Geosciences, University of Oslo
% Author: Alban Souche

% Computes the consistent tangent modulus for the Von Mises Elasto-plastic
% material with piece-wise linear isotropic hardening (plain strain)
% OLD VERSION BELOW derived from HYPLAS2 code (de Souza Neto, E. A.)


nel      = size(MESH.ELEMS,2);
D_el_ip  = zeros(nel,nip,9);


% Initialize entries to elastic tangent modulus
D_el_ip(:,:,1) = bsxfun(@times, 2*param.ep_G(MESH.elem_markers)+param.ep_lambda(MESH.elem_markers), ones(nip,1))' ;
D_el_ip(:,:,2) = bsxfun(@times, param.ep_lambda(MESH.elem_markers), ones(nip,1))';
D_el_ip(:,:,4) = bsxfun(@times, param.ep_lambda(MESH.elem_markers), ones(nip,1))';
D_el_ip(:,:,5) = bsxfun(@times, 2*param.ep_G(MESH.elem_markers)+param.ep_lambda(MESH.elem_markers), ones(nip,1))';
D_el_ip(:,:,9) = bsxfun(@times, param.ep_G(MESH.elem_markers), ones(nip,1))';



% Compute the elastiplastic tangent modulus at the itegration point in the plastic domain
if ~isequal(indx_p,[])
    
    % Ini. elastoplastic tangent matrix entries
    soid    = [1 1 0];
    foid    = [1 0 0;0 1 0; 0 0 .5];
    devprj  = foid - (soid'*soid)/3;
    soidmat = soid'*soid;
    R2G     = 2 * param.ep_G;
    R3G     = 3 * param.ep_G;
    Roo3D2  = sqrt(3/2);
    
    for ip = 1:nip
        
        indxel_ip = indx_p(ip,:);
        indxel_ip = find(indxel_ip); % collecting the relevant elements
        idel_ip   = MESH.elem_markers(indxel_ip);
        
        if indxel_ip
            
            S      = zeros(4,length(indxel_ip));
            S(1,:) = Sdev.xx(ip,indxel_ip);
            S(2,:) = Sdev.yy(ip,indxel_ip);
            S(3,:) = Sdev.xy(ip,indxel_ip);
            S(4,:) = Sdev.zz(ip,indxel_ip);
            
            d_g    = d_gamma(ip,indxel_ip);
            
            SNORM  = sqrt(S(1,:).^2 + S(2,:).^2 + S(4,:).^2 + 2*S(3,:).^2); % == sqrt(2*J2) by definition
            Q      = Roo3D2*SNORM;
            % Q = Q + R3G*d_gamma(ip,indxel_ip); % only use if the input deviatoric stresses are the updated ones but not sure about this implementation
            
            AFACT  = R2G(idel_ip).*(1-R3G(idel_ip).*d_g./Q);
            BFACT  = 6*(param.ep_G(idel_ip).^2) .* ( d_g./Q - 1./(R3G(idel_ip)+param.ep_H) ) ./ (SNORM.^2);
            
            Smat = [ S(1,:).*S(1,:); S(2,:).*S(1,:); S(3,:).*S(1,:); ...
                S(1,:).*S(2,:); S(2,:).*S(2,:); S(3,:).*S(2,:); ...
                S(1,:).*S(3,:); S(2,:).*S(3,:); S(3,:).*S(3,:)  ];  % Smat = S(1:3,1)*S(1:3,1)';  
            
            D_el_ip(indxel_ip,ip,:) = bsxfun(@plus, ( devprj(:)*AFACT + bsxfun(@times,Smat,BFACT) ), soidmat(:).*param.ep_K(idel_ip))';
            
        end
        
    end
    
end




% %%% OLD ELEMENT LOOP VERSION

% nel = size(MESH.ELEMS,2);
% D_el    = zeros(3, 3);
% D_el_ip3 = zeros(nel, nip, 9);
%
% % Ini. elastic tangent matrix entries
% D_el(1) = 2*param.mu+param.lambda;
% D_el(2) = param.lambda;
% D_el(4) = param.lambda;
% D_el(5) = 2*param.mu+param.lambda;
% D_el(9) = param.mu;
%
% % Ini. elastoplastic tangent matrix entries
% soid   = [1 1 0];
% foid   = [1 0 0;0 1 0; 0 0 .5];
% D_elpl = zeros(3, 3);
%
% % d_g = zeros(size(indx_p));
% % d_g(indx_p) = d_gamma;
% d_g = d_gamma;
%
% for iel = 1:nel
%
%     for ip = 1:nip
%
%         toeva = indx_p(ip,iel);
%
%         switch toeva
%
%             case 0 % elastic tangent matrix entries
%                 D_el_ip3(iel,ip,:) = D_el(:);
%
%             case 1 % elastoplastic tangent matrix entries
%
%                 devprj = zeros(3,3);
%                 for i = 1:3
%                     for j=1:3
%                         devprj(i,j) = foid(i,j) - ( soid(i)*soid(j) /3 );
%                     end
%                 end
%
%                 R2G    = 2 * param.G;
%                 R3G    = 3 * param.G;
%
%                 Roo3D2 = sqrt(3/2);
%
%                 S = zeros(4,1);
%                 S(1) = Sdev.xx(ip,iel);
%                 S(2) = Sdev.yy(ip,iel);
%                 S(3) = Sdev.xy(ip,iel);
%                 S(4) = Sdev.zz(ip,iel);
%                 SNORM  = sqrt(S(1)^2 + S(2)^2 + S(4)^2 + 2*S(3)^2);
%                 Q      = Roo3D2*SNORM;
%                 Qtrial = Q + R3G*d_g(ip,iel);
%
%                 AFACT  = R2G*(1-R3G*d_g(ip,iel)/Qtrial);
%                 BFACT  = 6*(param.G^2) * ( d_g(ip,iel)/Qtrial - 1/(R3G+param.H) ) / (SNORM^2);
%
%                 % upper triangle
%                 for i=1:3
%                     for j=i:3
%                         D_elpl(i,j) = AFACT*devprj(i,j) + BFACT*S(i)*S(j) + param.K*soid(i)*soid(j);
%                     end
%                 end
%
%                 % lower triangle
%                 for j=1:2
%                     for i=(j+1):3
%                         D_elpl(i,j) = D_elpl(j,i);
%                     end
%                 end
%
%
%                 D_el_ip3(iel,ip,:) = D_elpl(:);
%
%         end
%
%
%     end
% end

