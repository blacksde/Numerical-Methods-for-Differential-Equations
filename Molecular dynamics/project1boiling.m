%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% project 1 question 2
% boiling point

% set the initial value
Nx = 10;
Ny = 10;
Np = Nx*Ny;
Nt = 3000;  % Choose Nt = 3000

% set initial position

% for particles in odd rows, x = 0, 1, 2, 3...9
x1 = 0:(Nx-1);

% for particles in even rows, x = 0.5, 1.5, 2.5, 3.5 ...9.5
x2 = x1+0.5;

% for particle in odd rows, y = sqrt(3) * (0, 1, 2, 3, 4)
y1 = (0:(Ny/2-1))*sqrt(3);

% for particles in even rows, x = sqrt(3) * (0.5, 1.5, 2.5, 3.5. 4.5)
y2 = y1+sqrt(3)/2;

% coodinate for paricle in row 1,3,5,7,9
[X1,Y1] = meshgrid(x1,y1);
X1 = reshape(X1,1,Np/2);
Y1 = reshape(Y1,1,Np/2);

% coodinate for paricle in row 2,4,6,8,10
[X2,Y2] = meshgrid(x2,y2);
X2 = reshape(X2,1,Np/2);
Y2 = reshape(Y2,1,Np/2);

% combine them
x = [X1,X2];
y = [Y1,Y2];

% can use following line to plot the position
% plot(x,y,'o')

% decritize t
t = linspace(0,1,Nt);
h = t(2)-t(1);

% Test T from 0 to 3.
Tset = linspace(0,3,51);
D = zeros(1,length(Tset));

for n = 1:length(Tset)
    % show the term working on
    disp(num2str(n));
    
    % set initial u, v
    u = Tset(n)^0.5*randn(1,Np);
    v = Tset(n)^0.5*randn(1,Np);
    
    phye = [x;y;u;v];
    
    
    % solve equation using rk 4
    for j = 1:(length(t)-1)
        m1 = ffunction(phye);
        m2 = ffunction(phye+m1/2*h);
        m3 = ffunction(phye+m2/2*h);
        m4 = ffunction(phye+m3*h);
        phye = phye+h/6*(m1+2*m2+2*m3+m4);
    end
    
    % calculate distance
    for i = 1:Np
        % for each run calculate the distance with the ith particle
         x = phye(1,:);
         y = phye(2,:);
         
         % get coodinate of ith paricle
         xi = x(i);
         yi = y(i);
         
         % get coodinates of rest Np-1 paricles
         xo = x;
         xo(i)=[];
         yo = y;
         yo(i)=[];
         
         % distance
         r = sqrt((xo-xi).^2+(yo-yi).^2);
         D(n) = D(n)+sum(r);
    end
    
end

% divide the D by Np*(Np-1)
D = D/Np/(Np-1);

plot(Tset,D);
