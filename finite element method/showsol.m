function showsol(Elem,U)

%1D shape functions
N11=@(xi)0.5*(1-xi);
N12=@(xi)0.5*(1+xi);

%2D shape functions
N1=@(xi,eta) N11(xi).*N11(eta);
N2=@(xi,eta) N12(xi).*N11(eta);
N3=@(xi,eta) N12(xi).*N12(eta);
N4=@(xi,eta) N11(xi).*N12(eta);

hold on;
for ielem=1:length(Elem)
  xx=Elem(ielem).x;
  yy=Elem(ielem).y;
  GN=Elem(ielem).GN;

 [Xi, Eta]=meshgrid(-1:0.1:1, -1:0.1:1);
 
 X=N1(Xi,Eta)*xx(1)+N2(Xi,Eta)*xx(2)+N3(Xi,Eta)*xx(3)+N4(Xi,Eta)*xx(4);
 Y=N1(Xi,Eta)*yy(1)+N2(Xi,Eta)*yy(2)+N3(Xi,Eta)*yy(3)+N4(Xi,Eta)*yy(4);
 UU=N1(Xi,Eta)*U(GN(1))+N2(Xi,Eta)*U(GN(2))+N3(Xi,Eta)*U(GN(3))+N4(Xi,Eta)*U(GN(4));
 
 surf(X,Y,0*X,UU,'edgecolor','none');


 plot([xx xx(1)], [yy yy(1)],'k-');

 %show element numbers
 % text(mean(xx),mean(yy),num2str(ielem));

  for ii=1:4
      %show global numbers of the nodes
      % text(xx(ii),yy(ii),num2str(GN(ii)));
  end
  
end

colorbar;
hold off; 


end

