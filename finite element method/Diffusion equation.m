[Elem, Xnodes, Ynodes] = proj2mesh(6);

% get the number of nodes
Nnodes = length(Xnodes);
Nelem = length(Elem);

% construct K and M
K = zeros(Nnodes,Nnodes);
M = zeros(Nnodes,Nnodes);

for i = 1:Nelem
    
    X = Elem(i).x;
    Y = Elem(i).y;
    K(Elem(i).GN,Elem(i).GN) =K(Elem(i).GN,Elem(i).GN) + Kmatrix(X, Y);
    M(Elem(i).GN,Elem(i).GN) =M(Elem(i).GN,Elem(i).GN) + Mmatrix(X, Y);
    
end

% deal with the boundary condition
% find the location of left node 
locl = find(Xnodes == min(Xnodes));

M(locl,:) = 0;
K(locl,:) = 0;
M(locl,locl) = 1;

% find the location of right node 
locr = find(Xnodes == max(Xnodes));

M(locr,:) = 0;
K(locr,:) = 0;
M(locr,locr) = 1;

K1 = -M\K;

% now the question turn to the initial ODE
% U' = K1 * U
% U(left) = 1 & U(Domain\left) = 0

% set the initial condition

% set everywhere 0
U = zeros(Nnodes,1);
% set left 1
U(locl) = 1;

prob=2;%1 - backward Euler, 
       %2 - C-N,
       %3 - RK4
       
%discretize x coordinate
N=2e4;
x=linspace(0,5,N);

%step size
h=x(2)-x(1);

for n=1:N-1
    
  % backward Euler
   if prob==1
      %scheme
      U=(eye(Nnodes)-h*K1)\U;
      %multiplication factor
      G=@(z) 1./(1-z);
   end

   % C-N
   if prob==2
      %scheme
      U=(eye(Nnodes)-h/2*K1)\(eye(Nnodes)+h/2*K1)*U;    
      %multiplication factor
      G=@(z) (1+z/2)./(1-z/2);
   end
   
  % RK4
   if prob==3
     % scheme
      m1=K1*U;
      m2=K1*(U+m1/2*h);
      m3=K1*(U+m2/2*h);
      m4=K1*(U+m3*h);
      U=U+h/6*(m1+2*m2+2*m3+m4);
    %  multiplication factor
      G=@(z) 1+z+z.^2/2+z.^3/6+z.^4/24;
   end
   
   if(rem(n,100) == 0)
       disp(num2str(n));
   end
   
end

figure;
showsol(Elem,U);

