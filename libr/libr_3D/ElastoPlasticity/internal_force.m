function Fint = internal_force(MESH, INVJ, DetJ, Sdev, P, SOLVER) 

% Parameters
ndim         = size(MESH.NODES,1);
[nnodel,nel] = size(MESH.ELEMS);
nedof        = nnodel*ndim;
nip          = SOLVER.nip;
nelblo       = min(nel, SOLVER.nelblo);
nblo         = ceil(nel/nelblo);
il           = 1;
iu           = nelblo;

% Derivative of shape functions
[ipuv, ipw]  = ip_tetra(nip);
[ ~, dNdu]   = shp_deriv_tetra(ipuv, nnodel);

% Allocation
F_all   = zeros(nedof,nel); 
F_block = zeros(nelblo, nedof);

for ib = 1:nblo
    
    % Extract determinant of Jacobian
    detJ  = DetJ(il:iu);
    % Extract inverse of Jacobian
    invJx = INVJ(il:iu,1:3);
    invJy = INVJ(il:iu,4:6);
    invJz = INVJ(il:iu,7:9);
    % resert block vectors to zeros
    F_block(:) = 0;

    for ip=1:nip
        
        % Stresses at itegration points
        S1 = Sdev.xx(ip,il:iu)' + P(ip,il:iu)';
        S2 = Sdev.yy(ip,il:iu)' + P(ip,il:iu)';
        S3 = Sdev.zz(ip,il:iu)' + P(ip,il:iu)';
        S4 = Sdev.xy(ip,il:iu)';
        S5 = Sdev.yz(ip,il:iu)';
        S6 = Sdev.zx(ip,il:iu)';
        
        % Extract spatial derivative
        dNdui  = dNdu{ip};
        dNdx   = invJx*dNdui';
        dNdy   = invJy*dNdui';
        dNdz   = invJz*dNdui';
        
        % Weight of the itegration point
        weight = ipw(ip)*detJ;
       
        % Internal force vector
        F_block(:,1:ndim:nedof) =  F_block(:,1:ndim:nedof) + bsxfun(@times,(S1.*weight),dNdx) + bsxfun(@times,(S4.*weight),dNdy) + bsxfun(@times,(S6.*weight),dNdz);
        F_block(:,2:ndim:nedof) =  F_block(:,2:ndim:nedof) + bsxfun(@times,(S2.*weight),dNdy) + bsxfun(@times,(S4.*weight),dNdx) + bsxfun(@times,(S5.*weight),dNdz);
        F_block(:,3:ndim:nedof) =  F_block(:,3:ndim:nedof) + bsxfun(@times,(S3.*weight),dNdz) + bsxfun(@times,(S5.*weight),dNdy) + bsxfun(@times,(S6.*weight),dNdx);
        
    end
     
    F_all(:,il:iu)    = F_block';
    
    il  = il+nelblo;
    
    if(ib==nblo-1)
        nelblo 	= nel-iu;
        F_block = zeros(nelblo, nedof);
    end
    iu  = iu+nelblo;
end

ELEM_DOF = zeros(nedof, nel,'int32');
ELEM_DOF(1:ndim:end,:) = ndim*(MESH.ELEMS-1)+1;
ELEM_DOF(2:ndim:end,:) = ndim*(MESH.ELEMS-1)+2;
ELEM_DOF(3:ndim:end,:) = ndim*(MESH.ELEMS-1)+3;
Fint  = accumarray(ELEM_DOF(:), F_all(:));




%%% OLD ELEMENT LOOP VERSION
% 
% [ndim,nnod]=size(MESH.NODES);
% [nnodel,nel]=size(MESH.ELEMS);
% ndofel = nnodel*ndim;
% 
% nip          = size(S.xx,1);
% [ipuv, ipw] = ip_triangle(nip);
% [~,dNdU] = shape_functions_deriv('tri',nnodel,ipuv);
% 
% B=zeros((ndim+1)*ndim/2,ndofel);
% Fel=zeros(ndofel,nel);
% 
% for iel=1:nel
%     F = zeros(ndofel,1);
%     NODES_el=MESH.NODES(:,MESH.ELEMS(:,iel));
%     
%     for ip=1:nip
%           
%         Si = [S.xx(ip,iel); S.yy(ip,iel) ; S.xy(ip,iel)];
%         
%         dNdU_ip=dNdU(:,:,ip);
%         J=NODES_el*dNdU_ip;     %jacobian
%         det_J=det(J);           %determinant of the jacobian
%         invJ=inv(J);
%         dNdX=dNdU_ip*invJ;
%         for idim=1:ndim
%             B(idim,idim:ndim:end)=dNdX(:,idim)';
%         end
%         B(3,1:ndim:end)=dNdX(:,2)';
%         B(3,2:ndim:end)=dNdX(:,1)';
%             
%         F = F + ipw(ip)*det_J*(B'*Si);
%     end
%     Fel(:,iel) = F;
% end
% 
% tmp1 = ndim*(MESH.ELEMS-1)+1;
% tmp2 = ndim*(MESH.ELEMS-1)+2;
% tmp  = [tmp1(1,:); tmp2(1,:); tmp1(2,:); tmp2(2,:); tmp1(3,:); tmp2(3,:)];
% Fint = accumarray(tmp(:), Fel(:));





%    %%% LINEAR EXTRAPOLATION FROM EVALUATION PTS TO INTEGRATION PTS
% neva         = size(Sdev.xx,1);
% if neva==3
%     
%     [ieuv, ~]  = ip_triangle(neva);
%        
%     Ax = ieuv(1,1)*ones(1,nel);
%     Ay = ieuv(1,2)*ones(1,nel);
%     Bx = ieuv(2,1)*ones(1,nel);
%     By = ieuv(2,2)*ones(1,nel);
%     Cx = ieuv(3,1)*ones(1,nel);
%     Cy = ieuv(3,2)*ones(1,nel);
%     
%     pts_eva_loc = [0,0; 1,0; 0,1; .5,.5; 0,.5; .5,0; 1/3,1/3];
%     DDx = repmat(pts_eva_loc(:,1),1,nel);
%     DDy = repmat(pts_eva_loc(:,2),1,nel);
%     
%     Az = S.xx(1,:);
%     Bz = S.xx(2,:);
%     Cz = S.xx(3,:);
%     
%     a = (By-Ay).*(Cz-Az)-(Cy-Ay).*(Bz-Az);
%     b = (Bz-Az).*(Cx-Ax)-(Cz-Az).*(Bx-Ax);
%     c = (Bx-Ax).*(Cy-Ay)-(Cx-Ax).*(By-Ay);
%     d =-(a.*Ax + b.*Ay + c.*Az);
%     
%     Snodxx = -bsxfun(@rdivide, bsxfun( @plus,(bsxfun(@times,a,DDx) + bsxfun(@times,b,DDy)) , d) , c);
%     
%     Az = S.yy(1,:);
%     Bz = S.yy(2,:);
%     Cz = S.yy(3,:);
%     
%     a = (By-Ay).*(Cz-Az)-(Cy-Ay).*(Bz-Az);
%     b = (Bz-Az).*(Cx-Ax)-(Cz-Az).*(Bx-Ax);
%     c = (Bx-Ax).*(Cy-Ay)-(Cx-Ax).*(By-Ay);
%     d =-(a.*Ax + b.*Ay + c.*Az);
%     
%     Snodyy = -bsxfun(@rdivide, bsxfun( @plus,(bsxfun(@times,a,DDx) + bsxfun(@times,b,DDy)) , d) , c);
%     
%     Az = S.xy(1,:);
%     Bz = S.xy(2,:);
%     Cz = S.xy(3,:);
%     
%     a = (By-Ay).*(Cz-Az)-(Cy-Ay).*(Bz-Az);
%     b = (Bz-Az).*(Cx-Ax)-(Cz-Az).*(Bx-Ax);
%     c = (Bx-Ax).*(Cy-Ay)-(Cx-Ax).*(By-Ay);
%     d =-(a.*Ax + b.*Ay + c.*Az);
%     
%     Snodxy = -bsxfun(@rdivide, bsxfun( @plus,(bsxfun(@times,a,DDx) + bsxfun(@times,b,DDy)) , d) , c);
% 
% end



