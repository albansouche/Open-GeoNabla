function Disp = mechanical_solver(MESH, A_all, BC, rhs, SOLVER)

% Part of Open-GeoNabla, copyright GPLv3, 2018
% https://github.com/albansouche/Open-GeoNabla/
% Physics of Geological Processes (PGP) , The NJORD Centre, Dept of Geosciences, University of Oslo
% Author: Alban Souche


nnod = size(MESH.NODES,2);
sdof        = 2*nnod;

switch SOLVER.solver_type
    
    case 'direct_nonsym'
        
        opts_sparse.symmetric  = 0;
        opts_sparse.n_node_dof = 2;
        opts_sparse.nthreads   = SOLVER.nb_cpu;
        A   = sparse_create(MESH.ELEMS, A_all, opts_sparse); 
        clear A_all;
        % BC
        Disp         = zeros(sdof  , 1);
        Disp(BC.ind) = BC.val;
        Free         = 1:sdof;
        Free(BC.ind) = [];
        Rhs          = rhs - A(:,BC.ind)*Disp(BC.ind);
        A            = A(Free, Free);
        % Solving    
%         tic
        Disp(Free)   = A\Rhs(Free);
%         toc

    case 'direct_sym'
        opts_sparse.symmetric  = 1;
        opts_sparse.n_node_dof = 2;
        opts_sparse.nthreads   = SOLVER.nb_cpu;
        opts_sparse.cpu_affinity = 1;
        A = sparse_create(MESH.ELEMS, A_all, opts_sparse);
        clear A_all;
        % BC
        Disp         = zeros(sdof  , 1);
        Disp(BC.ind) = BC.val;
        Free         = 1:sdof;
        Free(BC.ind) = [];
        Rhs          = rhs - A*Disp - A'*Disp;
        A            = A(Free, Free);
        % Backslash
        A            = A + tril(A,-1)';
        Disp(Free)   = A\Rhs(Free);
        
        
    case 'direct_sym_chol'
        
        opts_sparse.symmetric  = 1;
        opts_sparse.n_node_dof = 2;
        opts_sparse.nthreads   = SOLVER.nb_cpu;
        opts_sparse.cpu_affinity = 1;
        A = sparse_create(MESH.ELEMS, A_all, opts_sparse);
        clear A_all;
        % BC
        Disp         = zeros(sdof  , 1);
        Disp(BC.ind) = BC.val;
        Free        = 1:sdof;
        Free(BC.ind)= [];
        Rhs         = rhs - A*Disp - A'*Disp;
        A           = A(Free, Free);

        % Factorizing and solving
%%%         [L,p] = lchol(A);
%%%%         F = factorize(A,'symmetric');
%%%%         Disp(Free) = F\Rhs(Free);

        [L,p] = chol(A, 'lower');
        if p>0
            fprintf(1,'Careful! Non positive definite striffness matrix. System solved with MATLAB backslash. \n');
            ttt
            A          = A + tril(A,-1)';
            Disp(Free) = A\Rhs(Free);
        else
            Disp(Free) = cs_ltsolve(L,cs_lsolve(L,Rhs(Free)));
        end
        
    case 'pcg'
% 
%         opts_sparse.symmetric  = 1;
%         opts_sparse.n_node_dof = 2;
%         opts_sparse.nthreads   = SOLVER.nb_cpu;
% %         opts_sparse.cpu_affinity = 1;
%         A = sparse_create(MESH.ELEMS, A_all, opts_sparse);
%         clear A_all;
%         % BC
%         Disp         = zeros(sdof  , 1);
%         Disp(BC.ind) = BC.val;
%         Free        = 1:sdof;
%         Free(BC.ind)= [];
%         Rhs         = rhs - A*Disp - A'*Disp;
%         Rhs         = Rhs(Free);
%         A           = A(Free, Free);
%          
% %         A = A + tril(A,-1)';
% %         [Disp(Free),~,~,iter] = pcg(A, Rhs(Free), 1e-9,10000);
% %         fprintf(1, ['PCG(',num2str(iter),'iter) ']);
% %         Disp(BC.ind) = BC.val;
% 
%         tic
%         precond = 1./spdiags(A,0);
%         A_conv = sparse_convert(A, opts_sparse);
%         Ax = @(x) spmv(A_conv,x);
%         [Disp(Free),~,~, iter] = pcg(Ax, Rhs, 1e-9, 10000, @(x) precond.*x);
%         fprintf(1, ['MPCG(',num2str(iter),'it) ']);toc
        

    case 'AGMG'
%         % AGMG
%         [Disp(Free),~,~,iter] = agmg(A,Rhs(Free),1,1e-9,100);
%         fprintf(1, ['AGMG(',num2str(iter),'it) ']);
%         Disp(BC.ind) = BC.val;

        
end





