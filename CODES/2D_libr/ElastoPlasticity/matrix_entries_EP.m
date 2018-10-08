function [A_all, Rhs_all, InvJ_all, DetJ_all] = matrix_entries_EP(MESH, D_el_ip, PARAM, SOLVER, elast)

% Part of Open-GeoNabla, copyright GPLv3, 2018
% https://github.com/albansouche/Open-GeoNabla/
% Physics of Geological Processes (PGP) , The NJORD Centre, Dept of Geosciences, University of Oslo
% Author: Alban Souche

%MODIFIED VERSION OF:
%MECHANICAL2D Two dimensional finite element mechanical problem solver of MILAMIN
%   Part of MILAMIN: MATLAB-based FEM solver for large problems
%   Version 1.0.1
%   Copyright (C) 2011, M. Dabrowski, M. Krotkiewski, D.W. Schmid
%   University of Oslo, Physics of Geological Processes
%   http://milamin.org
%   See License file for terms of use.

%==========================================================================
% MODEL INFO
%==========================================================================
[ndim,     ~] = size(MESH.NODES);
[nnodel, nel] = size(MESH.ELEMS);
nedof         = nnodel*ndim;
nip           = SOLVER.nip;

%==========================================================================
% BLOCKING PARAMETERS (nelblo must be < nel)
%==========================================================================
nelblo          = min(nel, SOLVER.nelblo);
nblo            = ceil(nel/nelblo);

%==========================================================================
% i) PREPARE INTEGRATION POINTS & DERIVATIVES wrt LOCAL COORDINATES
%==========================================================================
[IP_X, IP_w]    = ip_triangle(nip);
[   N, dNdu]    = shp_deriv_triangle(IP_X, nnodel);

%==================================================================
% DECLARE VARIABLES (ALLOCATE MEMORY)
%==================================================================
if or(elast==1, isequal(PARAM.ep_phi,PARAM.ep_xsi))
    A_block     = zeros(nelblo, nedof*(nedof+1)/2);
    A_all       = zeros(nedof*(nedof+1)/2,nel);
elseif max(PARAM.ep_phi)>min(PARAM.ep_xsi)
    A_block     = zeros(nelblo, nedof*nedof);
    A_all       = zeros(nedof*nedof,nel);
end
Rhs_all     = zeros(nedof,nel);
Rhs_block   = zeros(nelblo, nedof);
InvJ_all    = zeros(nel,4);
DetJ_all    = zeros(nel,1);

il          = 1;
iu          = nelblo;
%==================================================================
% i) BLOCK LOOP - MATRIX COMPUTATION
%==================================================================

for ib = 1:nblo
    %==============================================================
    % ii) FETCH DATA OF ELEMENTS IN BLOCK
    %==============================================================
    ECOORD_x = reshape( MESH.NODES(1,MESH.ELEMS(:,il:iu)), nnodel, nelblo);
    ECOORD_y = reshape( MESH.NODES(2,MESH.ELEMS(:,il:iu)), nnodel, nelblo);
    rho_el   = PARAM.ep_rho(MESH.elem_markers(il:iu))';
    
    %==============================================================
    % iii) INTEGRATION LOOP
    %==============================================================
    A_block(:)    = 0;
    Rhs_block(:)  = 0;
    

    for ip=1:nip

        %==========================================================
        % iv) LOAD SHAPE FUNCTIONS DERIVATIVES FOR INTEGRATION POINT
        %==========================================================
        Ni      =        N{ip};
        dNdui   =     dNdu{ip};

        %==========================================================
        % v) CALCULATE JACOBIAN, ITS DETERMINANT AND INVERSE
        %==========================================================
        Jx          = ECOORD_x'*dNdui;
        Jy          = ECOORD_y'*dNdui;
        detJ        = Jx(:,1).*Jy(:,2) - Jx(:,2).*Jy(:,1);
        
        invdetJ     = 1.0./detJ;
        invJx(:,1)  = +Jy(:,2).*invdetJ;
        invJx(:,2)  = -Jy(:,1).*invdetJ;
        invJy(:,1)  = -Jx(:,2).*invdetJ;
        invJy(:,2)  = +Jx(:,1).*invdetJ;
        
        %==========================================================
        % vi) DERIVATIVES wrt GLOBAL COORDINATES
        %==========================================================
        dNdx        = invJx*dNdui';
        dNdy        = invJy*dNdui';
        
        %==========================================================
        % vii) NUMERICAL INTEGRATION OF ELEMENT MATRICES
        %==========================================================
        weight      = IP_w(ip)*detJ;
                    
        % ------------------------A matrix-------------------------
        if elast==1 %elatic
            
            % Load elasto-plastic tangent matrix entries
            d11 = D_el_ip(il:iu,ip,1);
            d12 = D_el_ip(il:iu,ip,2);
            d33 = D_el_ip(il:iu,ip,9);
            
            indx  = 1;
            for i = 1:nnodel
                % x-velocity equation
                for j = i:nnodel
                    A_block(:,indx) = A_block(:,indx) + ( d11.*dNdx(:,i).*dNdx(:,j) + d33.*dNdy(:,i).*dNdy(:,j)).*weight;
                    indx = indx+1;
                    A_block(:,indx) = A_block(:,indx) + ( d12.*dNdx(:,i).*dNdy(:,j) + d33.*dNdy(:,i).*dNdx(:,j)).*weight;
                    indx = indx+1;
                end
                % y-velocity equation
                for j = i:nnodel
                    if(j>i)
                        A_block(:,indx) = A_block(:,indx) + ( d12.*dNdy(:,i).*dNdx(:,j) + d33.*dNdx(:,i).*dNdy(:,j)).*weight;
                        indx = indx+1;
                    end
                    A_block(:,indx) = A_block(:,indx) + ( d11.*dNdy(:,i).*dNdy(:,j) + d33.*dNdx(:,i).*dNdx(:,j)).*weight;
                    indx = indx+1;
                end
            end
            
        elseif PARAM.ep_phi==PARAM.ep_xsi % elastoplastic & associative plastic flow rule
            
            d11 = D_el_ip(il:iu,ip,1);
            d12 = D_el_ip(il:iu,ip,4);
            d13 = D_el_ip(il:iu,ip,7);
            d22 = D_el_ip(il:iu,ip,5);
            d23 = D_el_ip(il:iu,ip,8);
            d33 = D_el_ip(il:iu,ip,9);

            indx  = 1;
            for i = 1:nnodel
                
                % X-VELOCITY EQUATION
                for j = i:nnodel
                    A_block(:,indx) = A_block(:,indx) + ( d11.*dNdx(:,i).*dNdx(:,j) + d33.*dNdy(:,i).*dNdy(:,j) + ...
                        d13.*dNdx(:,i).*dNdy(:,j) + d13.*dNdx(:,j).*dNdy(:,i)).*weight;
                    indx = indx+1;
                    A_block(:,indx) = A_block(:,indx) + (d12.*dNdx(:,i).*dNdy(:,j) + d33.*dNdy(:,i).*dNdx(:,j) + ...
                        d13.*dNdx(:,i).*dNdx(:,j) + d23.*dNdy(:,j).*dNdy(:,i)).*weight;
                    indx = indx+1;
                end
                
                %Y-VELOCITY EQUATION
                for j = i:nnodel
                    if(j>i)
                        A_block(:,indx) = A_block(:,indx) + (d12.*dNdy(:,i).*dNdx(:,j) + d33.*dNdx(:,i).*dNdy(:,j) + ...
                            d13.*dNdx(:,i).*dNdx(:,j) + d23.*dNdy(:,j).*dNdy(:,i)).*weight;
                        indx = indx+1;
                    end
                    A_block(:,indx) = A_block(:,indx) + ( d22.*dNdy(:,i).*dNdy(:,j) + d33.*dNdx(:,i).*dNdx(:,j) + ...
                        d23.*dNdx(:,i).*dNdy(:,j) + d23.*dNdx(:,j).*dNdy(:,i)).*weight;
                    indx = indx+1;
                end
            end
            
        elseif PARAM.ep_phi>PARAM.ep_xsi  % elastoplastic & non-associative plastic flow rule
            
%             d11 = D_el_ip(il:iu,ip,1);
%             d12 = D_el_ip(il:iu,ip,4);
%             d13 = D_el_ip(il:iu,ip,7);
%             d21 = D_el_ip(il:iu,ip,2);
%             d22 = D_el_ip(il:iu,ip,5);
%             d23 = D_el_ip(il:iu,ip,8);
%             d31 = D_el_ip(il:iu,ip,3);
%             d32 = D_el_ip(il:iu,ip,6);
%             d33 = D_el_ip(il:iu,ip,9);
            
            d11 = D_el_ip(il:iu,ip,1);
            d12 = D_el_ip(il:iu,ip,2);
            d13 = D_el_ip(il:iu,ip,3);
            d21 = D_el_ip(il:iu,ip,4);
            d22 = D_el_ip(il:iu,ip,5);
            d23 = D_el_ip(il:iu,ip,6);
            d31 = D_el_ip(il:iu,ip,7);
            d32 = D_el_ip(il:iu,ip,8);
            d33 = D_el_ip(il:iu,ip,9);
            
            indx  = 1;
            for i = 1:nnodel
                % X-VELOCITY EQUATION
                for j = 1:nnodel
                    A_block(:,indx) = A_block(:,indx) + ( d11.*dNdx(:,i).*dNdx(:,j) + d31.*dNdx(:,j).*dNdy(:,i) + d13.*dNdx(:,i).*dNdy(:,j) + d33.*dNdy(:,i).*dNdy(:,j) ).*weight;
                    indx = indx+1;
                    A_block(:,indx) = A_block(:,indx) + ( d12.*dNdx(:,i).*dNdy(:,j) + d32.*dNdy(:,j).*dNdy(:,i) + d13.*dNdx(:,i).*dNdx(:,j) + d33.*dNdy(:,i).*dNdx(:,j) ).*weight;
                    indx = indx+1; 
                end
                for j = 1:nnodel
                %Y-VELOCITY EQUATION
                    A_block(:,indx) = A_block(:,indx) + ( d21.*dNdy(:,i).*dNdx(:,j) + d31.*dNdx(:,j).*dNdx(:,i) + d23.*dNdy(:,i).*dNdy(:,j) + d33.*dNdx(:,i).*dNdy(:,j) ).*weight;
                    indx = indx+1;
                    A_block(:,indx) = A_block(:,indx) + ( d22.*dNdy(:,i).*dNdy(:,j) + d32.*dNdy(:,j).*dNdx(:,i) + d23.*dNdy(:,i).*dNdx(:,j) + d33.*dNdx(:,i).*dNdx(:,j) ).*weight;
                    indx = indx+1;
                end
            end
            
        end
        
        % -----------------------Rhs vector------------------------
        Rhs_block(:,1:2:nedof) = Rhs_block(:,1:2:nedof) + PARAM.g(1)*(rho_el.*weight)*Ni';
        Rhs_block(:,2:2:nedof) = Rhs_block(:,2:2:nedof) + PARAM.g(2)*(rho_el.*weight)*Ni';
        
    end
    
    %==============================================================
    % ix) WRITE DATA INTO GLOBAL STORAGE
    %==============================================================
    A_all(:,il:iu)    = A_block';
    Rhs_all(:,il:iu)  = Rhs_block';
    InvJ_all(il:iu,:) = [invJx invJy];
    DetJ_all(il:iu,:) = detJ;
    
    %==============================================================
    % READJUST START, END AND SIZE OF BLOCK. REALLOCATE MEMORY
    %==============================================================
    il  = il+nelblo;
    if(ib==nblo-1)
        nelblo      = nel-iu;
        if or(elast==1, isequal(PARAM.ep_phi,PARAM.ep_xsi))
            A_block     = zeros(nelblo, nedof*(nedof+1)/2);
        elseif max(PARAM.ep_phi)>min(PARAM.ep_xsi)
            A_block     = zeros(nelblo, nedof*nedof);
        end
        Rhs_block   = zeros(nelblo, nedof);
        invJx       = zeros(nelblo, ndim);
        invJy       = zeros(nelblo, ndim);
        detJ        = zeros(nelblo, 1);
    end
    iu  = iu+nelblo;
end
