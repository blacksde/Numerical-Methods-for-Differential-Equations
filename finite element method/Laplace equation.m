[Elem, Xnodes, Ynodes] = proj2mesh(6);

% get the number of nodes
Nnodes = length(Xnodes);
Nelem = length(Elem);

% construct K and F
K = zeros(Nnodes,Nnodes);
F = zeros(Nnodes,1);

for i = 1:Nelem
    X = Elem(i).x;
    Y = Elem(i).y;
    K(Elem(i).GN,Elem(i).GN) =K(Elem(i).GN,Elem(i).GN) + Kmatrix(X, Y);
end

% deal with the boundary condition
% find the location of left node 
locl = find(Xnodes == min(Xnodes));

K(locl,:) = 0;
K(locl,locl) = 1;
F(locl) = 1;

% find the location of right node 
locr = find(Xnodes == max(Xnodes));

K(locr,:) = 0;
K(locr,locr) = 1;
F(locr) = 0;

% solve U by inv(K)*F
U = K\F;

figure;
showsol(Elem,U);

