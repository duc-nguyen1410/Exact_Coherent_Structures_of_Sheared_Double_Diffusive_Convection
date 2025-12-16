function [eigenvector,eigenvalue] = eig_DDC(Ra,Pr,Lambda,tau,kx,ky,N)
     % This function computes eigenvalues of linearized SDDC equations
     
     I=eye(N,N);
     O=zeros(N,N);

     % Differentiation Matrices 
     [z, D] = chebdif(N, 2);
     D1 = D(:,:,1); 
     D2 = D(:,:,2);

     % interval transformation [-1,1] to [-0.5,0.5]
     z = z/2;
     D1 = 2*D1;
     D2 = 4*D2;

     Ub = 0 * z;

     DX = 1i*kx*I;
     DY = 1i*ky*I;
     DZ = D1;
     Laplacian = -(kx^2+ky^2)*I + D2;

     Mu = sqrt(Pr/Ra)*Laplacian;
     Mv = sqrt(Pr/Ra)*Laplacian;
     Mw = sqrt(Pr/Ra)*Laplacian;
     Mt = (1.0/sqrt(Pr*Ra))*Laplacian;
     Mc = (tau/sqrt(Pr*Ra))*Laplacian;
     G = diag(-Lambda)*I;


     A = [Mu, O, O,-DX, O, O;%u
           O,Mv, O,-DY, O, O;%v
           O, O,Mw,-DZ, I, G;%w
          DX,DY,DZ,  O, O, O;%p
           O, O, I,  O,Mt, O;%t
           O, O, I,  O, O,Mc];%c
     B = [I, O, O, O, O, O;%u
          O, I, O, O, O, O;%v
          O, O, I, O, O, O;%w
          O, O, O, O, O, O;%p
          O, O, O, O, I, O;%t
          O, O, O, O, O, I];%c

     % Apply Dirichlet BCs at top/bottom for u, v, w, t, c (but not p)
     for var = 0:5
          if var == 3  % skip p
               continue
          end
          row_top = 1 + var * N;
          row_bot = N + var * N;
          A(row_top,:) = 0; A(row_bot,:) = 0;
          B(row_top,:) = 0; B(row_bot,:) = 0;
          A(row_top, row_top) = 1;
          A(row_bot, row_bot) = 1;
          %    B(row_top, row_top) = 0;
          %    B(row_bot, row_bot) = 0;
     end

     [eigenvector,eigenvalue] = eig(A,B);

     % Optional: sort eigenvalues by real part
     % sigma = diag(eigenvalue);
     % [~, idx] = sort(real(sigma), 'descend');
     % eigenvalue = diag(sigma(idx));
     % eigenvector = eigenvector(:, idx);
end