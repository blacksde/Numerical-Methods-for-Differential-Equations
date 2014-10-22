%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this is the function for the Mmatrix of every element 
% w.r.t to Parabolic and Hyperbolic function
% X is the x cood of node with length 4
% Y is the y cood of node with length 4
% M is a Mmatrix w.r.t the element defined by X and Y


function  M= Mmatrix(X, Y)

% to check the lengthes of X and Y are 4
if(length(X) ~= 4 || length(Y) ~= 4)
     error('dim unfit');
end

%Gauss points
xgp=[-sqrt(3/5), 0, sqrt(3/5)];

%weights
wg=[5/9, 8/9, 5/9];

% introduce N vector of Isoparameteric mapping
N = @(Xi,Eta)(0.25*[(1-Xi)*(1-Eta),(1+Xi)*(1-Eta),(1+Xi)*(1+Eta),(1-Xi)*(1+Eta)]);

% introduce J matrix of Isoparameteric mapping
J = @(xi,eta)(0.25*[-1+eta,1-eta,1+eta,-1-eta; -1+xi, -1-xi,1+xi,1-xi]*transpose([X;Y]));

M = zeros(4,4);

for i = 1:length(xgp)
    for j = 1:length(xgp)
       Jnode = J(xgp(i),xgp(j));
       Nnode = N(xgp(i),xgp(j));         
        
       M = M +(abs(det(Jnode))* wg(i)*wg(j))*transpose(Nnode)*Nnode;
       
    end
end

end