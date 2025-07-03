% function [U, A] = computePDEModel(U, A, Ic, A_Ic,  dc1_dx, dc2_dx, dc1_dy, dc2_dy,params)
function [U, A, U_a, A_a] = computePDEModel_op(U, A, U_a, A_a, Ic, Inds_Ic, A_Ic,  dc1_dx, dc2_dx, dc1_dy, dc2_dy,params)
    % Copy the current state of u
    % tic;
    U_old = U;   % U_old(Ic == 1) = 0;
    A_old = A;   % A_old(Ic == 1) = 0;
    U_X_old = U_old;
    U_Y_old = U_old;
    A_X_old = A_old;
    A_Y_old = A_old;
    
    A_a_old = A_a;% A_a_old(Ic == 1) = 0;
    A_a_X_old = A_a_old;
    A_a_Y_old = A_a_old;
    % disp(['Initial store: ', num2str(toc), ' seconds.'])
    % tic;
    % for ii=2:params.Ny-1
    %     for jj=2:params.Nx-1
    %         if Ic(ii,jj) == 1
    %             if Ic(ii+1,jj) == 0
    %                 U_X_old(ii,jj) =      U(ii+1,jj);
    %                 A_X_old(ii,jj) =      A(ii+1,jj);
    %                 A_a_X_old(ii,jj) =  A_a(ii+1,jj);
    %             % end
    %             elseif Ic(ii-1,jj) == 0
    %                 U_X_old(ii,jj) =      U(ii-1,jj);
    %                 A_X_old(ii,jj) =      A(ii-1,jj);
    %                 A_a_X_old(ii,jj) =  A_a(ii-1,jj);
    %             end
    % 
    %             if Ic(ii,jj+1) == 0
    %                 U_Y_old(ii,jj) =      U(ii,jj+1);
    %                 A_Y_old(ii,jj) =      A(ii,jj+1);
    %                 A_a_Y_old(ii,jj) =  A_a(ii,jj+1);
    %             % end
    %             elseif Ic(ii,jj-1) == 0
    %                 U_Y_old(ii,jj) =      U(ii,jj-1);
    %                 A_Y_old(ii,jj) =      A(ii,jj-1);
    %                 A_a_Y_old(ii,jj) =  A_a(ii,jj-1);
    %             end
    %         end
    %     end
    % end
    % for kk = 1:length(Inds_Ic)
    %     ii = Inds_Ic(kk, 1);
    %     jj = Inds_Ic(kk, 2);
    %     if Ic(ii+1,jj) == 0
    %         U_X_old(ii,jj) =      U(ii+1,jj);
    %         A_X_old(ii,jj) =      A(ii+1,jj);
    %         A_a_X_old(ii,jj) =  A_a(ii+1,jj);
    %     % end
    %     elseif Ic(ii-1,jj) == 0
    %         U_X_old(ii,jj) =      U(ii-1,jj);
    %         A_X_old(ii,jj) =      A(ii-1,jj);
    %         A_a_X_old(ii,jj) =  A_a(ii-1,jj);
    %     end
    % 
    %     if Ic(ii,jj+1) == 0
    %         U_Y_old(ii,jj) =      U(ii,jj+1);
    %         A_Y_old(ii,jj) =      A(ii,jj+1);
    %         A_a_Y_old(ii,jj) =  A_a(ii,jj+1);
    %     % end
    %     elseif Ic(ii,jj-1) == 0
    %         U_Y_old(ii,jj) =      U(ii,jj-1);
    %         A_Y_old(ii,jj) =      A(ii,jj-1);
    %         A_a_Y_old(ii,jj) =  A_a(ii,jj-1);
    %     end
    % end
    Ic_iip1_jj = Inds_Ic.Ic_iip1_jj;
    iip1_jj = Inds_Ic.iip1_jj;

    Ic_iim1_jj = Inds_Ic.Ic_iim1_jj;
    iim1_jj = Inds_Ic.iim1_jj;

    Ic_ii_jjp1 = Inds_Ic.Ic_ii_jjp1;
    ii_jjp1 = Inds_Ic.ii_jjp1;

    Ic_ii_jjm1 = Inds_Ic.Ic_ii_jjm1;
    ii_jjm1 = Inds_Ic.ii_jjm1;
    
    U_X_old(iip1_jj) =    U(Ic_iip1_jj);
    A_X_old(iip1_jj) =    A(Ic_iip1_jj);
    A_a_X_old(iip1_jj) =  A_a(Ic_iip1_jj);

    U_X_old(iim1_jj) =    U(Ic_iim1_jj);
    A_X_old(iim1_jj) =    A(Ic_iim1_jj);
    A_a_X_old(iim1_jj) =  A_a(Ic_iim1_jj);

    U_Y_old(ii_jjp1) =    U(Ic_ii_jjp1);
    A_Y_old(ii_jjp1) =    A(Ic_ii_jjp1);
    A_a_Y_old(ii_jjp1) =  A_a(Ic_ii_jjp1);

    U_Y_old(ii_jjm1) =    U(Ic_ii_jjm1);
    A_Y_old(ii_jjm1) =    A(Ic_ii_jjm1);
    A_a_Y_old(ii_jjm1) =  A_a(Ic_ii_jjm1);

    % disp(['1st boundary conditions: ', num2str(toc), ' seconds.'])
    % tic;

    

    % Reflective edges
    % u_X_extended = padarray(U_X_old, [1, 1], 'replicate'); % Extend boundaries by mirroring
    % u_Y_extended = padarray(U_Y_old, [1, 1], 'replicate'); % Extend boundaries by mirroring
    u_X_extended = myPadArray(U_X_old);
    u_Y_extended = myPadArray(U_Y_old);
    % disp(['U array padding: ', num2str(toc), ' seconds.'])

    % pad_size = length(U_X_old)+2;
    % pad1 = 3:pad_size;
    % pad2 = 2:(pad_size-1);
    % pad3 = 1:(pad_size-2);
    % Extract shifted matrices (inner region only)
    u_shift_xp = u_X_extended(3:end, 2:end-1); % Shift down
    u_shift_xm = u_X_extended(1:end-2, 2:end-1); % Shift up
    u_shift_yp = u_Y_extended(2:end-1, 3:end); % Shift right
    u_shift_ym = u_Y_extended(2:end-1, 1:end-2); % Shift left
    % u_shift_xp = u_X_extended(pad1, pad2); % Shift down
    % u_shift_xm = u_X_extended(pad3, pad2); % Shift up
    % u_shift_yp = u_Y_extended(pad2, pad1); % Shift right
    % u_shift_ym = u_Y_extended(pad2, pad3); % Shift left
    % disp(['U slicing: ', num2str(toc), ' seconds.'])

    % Central difference for diffusion
    uxx = (u_shift_xp - 2 * U_X_old + u_shift_xm) / (params.dx)^2;
    uyy = (u_shift_yp - 2 * U_Y_old + u_shift_ym) / (params.dy)^2;
    diffusion = params.D * (uxx + uyy);
    
    % For du_dx and du_dy terms
    du1_dx = (U_X_old + u_shift_xp) / 2;
    du2_dx = (U_X_old + u_shift_xm) / 2;
    du1_dy = (U_Y_old + u_shift_yp) / 2;
    du2_dy = (U_Y_old + u_shift_ym) / 2;
    % Taxis terms
    taxis_x = (du1_dx .* dc1_dx - du2_dx .* dc2_dx) / params.dx;
    taxis_y = (du1_dy .* dc1_dy - du2_dy .* dc2_dy) / params.dy;
    taxis = params.C_chi * (taxis_x + taxis_y);
    
    Activation_function = 1./(1+exp(-100.* ( A_old - params.activatedAge )));
    % Activation_function = 1./(1+exp(-100.* ( A_old - (params.activatedAge).*U_old )));
    U = U_old + params.dt * (diffusion - taxis);% - params.dt * U_old.*(Activation_function ) ;

    % U_a = U_a + params.dt*U_old.*(Activation_function );
    
    % Update antigen accumulation
    % Reflective edges
    % A_X_extended = padarray((A_X_old), [1, 1], 'replicate'); % Extend boundaries by mirroring
    % A_Y_extended = padarray((A_Y_old), [1, 1], 'replicate'); % Extend boundaries by mirroring
    A_X_extended = myPadArray(A_X_old);
    A_Y_extended = myPadArray(A_Y_old);

    % Extract shifted matrices (inner region only)
    A_shift_xp = A_X_extended(3:end, 2:end-1); % Shift down
    A_shift_xm = A_X_extended(1:end-2, 2:end-1); % Shift up
    A_shift_yp = A_Y_extended(2:end-1, 3:end); % Shift right
    A_shift_ym = A_Y_extended(2:end-1, 1:end-2); % Shift left
    % A_shift_xp = A_X_extended(pad1, pad2); % Shift down
    % A_shift_xm = A_X_extended(pad3, pad2); % Shift up
    % A_shift_yp = A_Y_extended(pad2, pad1); % Shift right
    % A_shift_ym = A_Y_extended(pad2, pad3); % Shift left
    
    % Central difference for diffusion
    Axx = (A_shift_xp - 2 * (A_X_old) + A_shift_xm) / (params.dx)^2;
    Ayy = (A_shift_yp - 2 * (A_Y_old) + A_shift_ym) / (params.dy)^2;
    A_diffusion = params.D * (Axx + Ayy);
    
    % For du_dx and du_dy terms
    dA1_dx = ((A_X_old) + A_shift_xp) / 2;
    dA2_dx = ((A_X_old) + A_shift_xm) / 2;
    dA1_dy = ((A_Y_old) + A_shift_yp) / 2;
    dA2_dy = ((A_Y_old) + A_shift_ym) / 2;
    % Taxis terms
    taxis_x = (dA1_dx .* dc1_dx - dA2_dx .* dc2_dx) / params.dx;
    taxis_y = (dA1_dy .* dc1_dy - dA2_dy .* dc2_dy) / params.dy;
    A_taxis = params.C_chi * (taxis_x + taxis_y);
    

    A = A_old + params.dt * ( A_diffusion - A_taxis  + params.P_A .* A_Ic .* (1./(params.activatedAge)).*(U_old - (A_old+A_a_old) )  ...
                       - params.P_D .* (1.0 - A_Ic) .* A_old/params.activatedAge ) ...
                       - params.dt .* params.P_A .* (1./(params.activatedAge^2)) .* A_old.*A_Ic;

    % Update antigen accumulation
    % Reflective edges
    % A_X_extended = padarray((A_a_X_old), [1, 1], 'replicate'); % Extend boundaries by mirroring
    % A_Y_extended = padarray((A_a_Y_old), [1, 1], 'replicate'); % Extend boundaries by mirroring
    A_X_extended = myPadArray(A_a_X_old);
    A_Y_extended = myPadArray(A_a_Y_old);
    % Extract shifted matrices (inner region only)
    A_shift_xp = A_X_extended(3:end, 2:end-1); % Shift down
    A_shift_xm = A_X_extended(1:end-2, 2:end-1); % Shift up
    A_shift_yp = A_Y_extended(2:end-1, 3:end); % Shift right
    A_shift_ym = A_Y_extended(2:end-1, 1:end-2); % Shift left
    % A_shift_xp = A_X_extended(pad1, pad2); % Shift down
    % A_shift_xm = A_X_extended(pad3, pad3); % Shift up
    % A_shift_yp = A_Y_extended(pad2, pad1); % Shift right
    % A_shift_ym = A_Y_extended(pad2, pad3); % Shift left
    
    % Central difference for diffusion
    Axx = (A_shift_xp - 2 * (A_a_X_old) + A_shift_xm) / (params.dx)^2;
    Ayy = (A_shift_yp - 2 * (A_a_Y_old) + A_shift_ym) / (params.dy)^2;
    A_diffusion = params.D * (Axx + Ayy);
    
    % For du_dx and du_dy terms
    dA1_dx = ((A_a_X_old) + A_shift_xp) / 2;
    dA2_dx = ((A_a_X_old) + A_shift_xm) / 2;
    dA1_dy = ((A_a_Y_old) + A_shift_yp) / 2;
    dA2_dy = ((A_a_Y_old) + A_shift_ym) / 2;
    % Taxis terms
    taxis_x = (dA1_dx .* dc1_dx - dA2_dx .* dc2_dx) / params.dx;
    taxis_y = (dA1_dy .* dc1_dy - dA2_dy .* dc2_dy) / params.dy;
    A_taxis = params.C_chi * (taxis_x + taxis_y);
    
    A_a = A_a_old + params.dt .* params.P_A .* (1./(params.activatedAge^2)) .* A_old.*A_Ic + params.dt * ( A_diffusion - A_taxis) ;
    % disp(['Matrix solutions: ', num2str(toc), ' seconds.'])
    % tic;
    Ic_vals = Inds_Ic.Ic_vals;
    A(Ic_vals) = 0;
    U(Ic_vals) = 0;
    A_a(Ic_vals) = 0;
    % disp(['2nd boundary conditions: ', num2str(toc), ' seconds.'])

end

%%
function paddedM = myPadArray(M)
    paddedM = zeros(size(M)+2);
    
    paddedM(1,2:(end-1)) = M(1,:);
    paddedM(end,2:(end-1)) = M(end,:);
    
    paddedM(2:(end-1),1) = M(:,1);
    paddedM(2:(end-1),end) = M(:,end);
    
    paddedM(2:(end-1),2:(end-1)) = M(:,:);
    
end