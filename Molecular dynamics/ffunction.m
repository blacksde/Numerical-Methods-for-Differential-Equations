%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function given by (5) on the project page
% phy stand for [x,y,u,v]
% return dphy is the [dx,dy,du,dv]


function dphy = ffunction(phy)

[Nphy, Np] = size(phy);

% if the row number is not 4,return error
if(Nphy ~= 4)
    error('dim unfit');
end

x = phy(1,:);
y = phy(2,:);
u = phy(3,:);
v = phy(4,:);
% seperate phy to x, y, u, v

% dx is u and dy is v
dx = u;
dy = v;

% construct vector du, dv
du = zeros(1,Np);
dv = zeros(1,Np);

for i = 1:Np
    % for each run calculate the ith particle
    
    % get coodinate of ith paricle
    xi = x(i);
    yi = y(i);
    
    % get coodinates of rest Np-1 paricles
    xo = x;
    xo(i)=[];
    yo = y;
    yo(i)=[];
    
    % distancs between them
    r = sqrt((xo-xi).^2+(yo-yi).^2);
    
    % du = -sum(dphi(r_ij)*(x_i-x_j)/r_ij)
    du(i) = -sum(dpotential(r).*(xi-xo)./r);
    dv(i) = -sum(dpotential(r).*(yi-yo)./r);
end

% return dphy as [dx,dy,du,dv]
dphy = [dx;dy;du;dv];

end