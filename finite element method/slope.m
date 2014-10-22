% MM = 9;
% this is for diffusion equation

MM = 7;
% this is for wave equation
% less than above to pretend out of memory

lam = zeros(1,MM);


for n = 1:MM
[Elem, Xnodes, Ynodes] = proj2mesh(n);

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

% K0 for diffusion equation
K0 = -M\K;

% K1 for wave equation
K1 = [zeros(Nnodes), eye(Nnodes);K0,zeros(Nnodes)];

%lam(n) = max(abs(eig(K0)));

lam(n) = max(abs(eig(K1)));



end

figure;
plot(-log(1:MM), -log(lam));

[slope, interrept] = leastsquare(-log(1:MM), -log(lam));