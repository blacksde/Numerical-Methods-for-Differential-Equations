%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% MATH 405 PROJECT 1 %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

% set initial condition of u and v
u = T^0.5*randn(1,Np);
v = T^0.5*randn(1,Np);

phy = [x;y;u;v];

% solve the "real" sulution
Nt = 5000;
phyr = phy;
t = linspace(0,5,Nt);
h = t(2)-t(1);

for j = 1:(length(t)-1)
    % rk 4 
    m1 = ffunction(phyr);
    m2 = ffunction(phyr+m1/2*h);
    m3 = ffunction(phyr+m2/2*h);
    m4 = ffunction(phyr+m3*h);
        
    phyr = phyr+h/6*(m1+2*m2+2*m3+m4);
    disp(num2str(j));
    
end

% solve the estimate sulution
Nt =  round(1000*1.1.^(1:15));
error = zeros(1,length(Nt));
for i = 1:length(Nt)
    
    % rk 4 
    phye = phy;
    t = linspace(0,5,Nt(i));
    h = t(2)-t(1);
    
    for j = 1:(length(t)-1)
        m1 = ffunction(phye);
        m2 = ffunction(phye+m1/2*h);
        m3 = ffunction(phye+m2/2*h);
        m4 = ffunction(phye+m3*h);
        phye = phye+h/6*(m1+2*m2+2*m3+m4);
    end
    
    % get the error with the exact solution
    error(i) = sqrt(sum(sum((phye-phyr).^2))/Np);
    disp(num2str(Nt(i)));
    
end

figure;
plot(log(Nt), log(error));