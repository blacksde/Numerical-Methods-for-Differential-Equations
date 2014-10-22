[Elem, Xnodes, Ynodes] = proj2mesh(6);

prob = 3; % 2 for diffusion equation
          % 3 for wave equation
method = 2; % 1 for Backward Euler
            % 2 for C-N
            % 3 for RK 4

if prob == 2
    N = [1,2,4,10]*5e3;
end

if prob == 3
    N = [1,2,4,10]*1e2;
end

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

if prob == 2
% for diffusion equation
K1 = -M\K;
end

if prob == 3
% for wave equation
K0 = -M\K;
K1 = [zeros(Nnodes), eye(Nnodes);K0,zeros(Nnodes)];
end


%  stability analysis
xmax=6;
i=sqrt(-1);
[xx,yy]=meshgrid(-xmax:0.01:xmax,-xmax:0.01:xmax);


if method == 1
    % backward Euler
    G=@(z) 1./(1-z);
elseif method == 2
        % C-N
    G=@(z) (1+z/2)./(1-z/2);
elseif method == 3
    % RK 4
    G=@(z) 1+z+z.^2/2+z.^3/6+z.^4/24;
end

GG=abs(G(xx+i*yy));
GG(GG>1)=NaN;

figure;
for j = 1:4
%discretize x coordinate
x=linspace(0,5,N(j));

%step size
h=x(2)-x(1);
%plot the stability region
subplot(2,2,j);
contourf(xx,yy,GG,50,'edgecolor','none');
hold on;

%calculate the eigenvalues
lambda=eig(K1);

%show the eigenvalues
plot(real(lambda)*h,imag(lambda)*h,'ko','markerfacecolor','k');

axis equal;
axis([-xmax xmax -xmax xmax]);

title(['N=',num2str(N(j))], 'FontWeight','bold')
xlabel('Re(\lambda h)','fontsize',12);
ylabel('Im(\lambda h)','fontsize',12);
disp(num2str(j));
end