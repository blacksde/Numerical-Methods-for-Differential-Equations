%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this is the function for the Kmatrix of every element w.r.t to 
% seconde order partial differential function
% X is the x cood of node with length 4
% Y is the y cood of node with length 4
% K is a matrix w.r.t the element defined by X and Y

function  K= Kmatrix(X, Y)


% to check the lengthes of X and Y are 4
if(length(X) ~= 4 || length(Y) ~= 4)
     error('dim unfit');
end

%Gauss points
xgp=[-sqrt(3/5), 0, sqrt(3/5)];

%weights
wg=[5/9, 8/9, 5/9];

% introduce B matrix base on domain after Iso mapping
B = @(xi,eta)(0.25*[-1+eta,1-eta,1+eta,-1-eta; -1+xi, -1-xi,1+xi,1-xi]);

% introduce J matrix of Isoparameteric mapping
J = @(xi,eta)(0.25*[-1+eta,1-eta,1+eta,-1-eta; -1+xi, -1-xi,1+xi,1-xi]*transpose([X;Y]));

% set initial value of K

K = zeros(4,4);

% calculate intergrals by 2D G-L method with node xgp and weight wg

for i = 1:length(xgp)
    for j = 1:length(xgp)
       Jnode = J(xgp(i),xgp(j));
       Bnode = B(xgp(i),xgp(j));         
        
       K = K +(abs(det(Jnode))* wg(i)*wg(j))*transpose(Jnode\Bnode)*(Jnode\Bnode);
    end
end

end