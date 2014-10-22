function [Elem, Xnodes, Ynodes] = proj2mesh(Nl)

Nl=Nl+1;
Nphi=Nl*8+1;

xmult=3;

%define outer and inner layers
t1=linspace(0,pi,Nphi);
t2=linspace(pi/16,15/16*pi,Nphi);

aa=sin(t2(1));

xoutr=xmult*(sin(t2)-aa);
youtr=cos(2*t2+pi/2);
xintr=xmult*(0.3+0.5*sin(t1)-aa);
yintr=0.5*cos(2*t1+pi/2);

xoutl=xmult*(sin(t2+pi)+aa);
youtl=cos(2*t2+pi/2+pi);
xintl=xmult*(-0.3+0.5*sin(t1+pi)+aa);
yintl=0.5*cos(2*t1+pi/2+pi);


%coordinates of the layers
Xl=zeros(Nl,Nphi-2);
Yl=zeros(Nl,Nphi-2);
Xr=zeros(Nl,Nphi-2);
Yr=zeros(Nl,Nphi-2);

%coordinates of the center piece
xcoord=[(-0.3+aa)*xmult 0 (0.3-aa)*xmult 0];
ycoord=[0 -cos(2*pi/16+pi/2) 0 cos(2*pi/16+pi/2)];
Xc=zeros(Nl,Nl);
Yc=zeros(Nl,Nl);
ind=2:length(xoutr)-1;

for in=1:Nl
   xr=xoutr(ind)*(1-(in-1)/(Nl-1))+xintr(ind)*(in-1)/(Nl-1);
   yr=youtr(ind)*(1-(in-1)/(Nl-1))+yintr(ind)*(in-1)/(Nl-1);
   xl=xoutl(ind)*(1-(in-1)/(Nl-1))+xintl(ind)*(in-1)/(Nl-1);
   yl=youtl(ind)*(1-(in-1)/(Nl-1))+yintl(ind)*(in-1)/(Nl-1);
  
   Xl(in,:)=xl;   
   Yl(in,:)=yl; 
   Xr(in,:)=xr;   
   Yr(in,:)=yr; 
   
   xx=linspace(xcoord(2), xcoord(1),Nl);
   dxx=linspace(xcoord(3)-xcoord(2), xcoord(4)-xcoord(1),Nl)/(Nl-1);
   yy=linspace(ycoord(2), ycoord(1),Nl);
   dyy=linspace(ycoord(3)-ycoord(2), ycoord(4)-ycoord(1),Nl)/(Nl-1);
  
   Xc(in,:)=xx+(in-1)*dxx;   
   Yc(in,:)=yy+(in-1)*dyy;
end



%number of nodes
Nnodes=Nl*(2*Nphi-4+Nl);
%number of elements
Nelem=(Nl-1)*(2*Nphi-3+Nl);


%location of the nodes
Xnodes=zeros(1,Nnodes);
Ynodes=zeros(1,Nnodes);

%properties of all elements
Elem(Nelem)=struct('x',[],'y',[],'GN',[]);


%left loop
for in =1:Nl-1
    
    for iphi=1:Nphi-3
        ielem=iphi+(in-1)*(Nphi-3);
        
        GN1=(in-1)*(Nphi-2)+iphi;
        GN2=(in-1)*(Nphi-2)+iphi+1;
        GN3=(in)*(Nphi-2)+iphi+1;
        GN4=(in)*(Nphi-2)+iphi;

        Xnodes([GN1 GN2 GN3 GN4])=[Xl(in,iphi) Xl(in,iphi+1) Xl(in+1,iphi+1) Xl(in+1,iphi)];
        Ynodes([GN1 GN2 GN3 GN4])=[Yl(in,iphi) Yl(in,iphi+1) Yl(in+1,iphi+1) Yl(in+1,iphi)];
        Elem(ielem).x=Xnodes([GN1 GN2 GN3 GN4]);
        Elem(ielem).y=Ynodes([GN1 GN2 GN3 GN4]);
        Elem(ielem).GN=[GN1 GN2 GN3 GN4];
    end

end
Nelemleft=ielem;
GNleft=Nl*(Nphi-2);

%right loop
for in =1:Nl-1
    
    for iphi=1:Nphi-3
        ielem=Nelemleft+iphi+(in-1)*(Nphi-3);
        
        GN1=GNleft+(in-1)*(Nphi-2)+iphi;
        GN2=GNleft+(in-1)*(Nphi-2)+iphi+1;
        GN3=GNleft+(in)*(Nphi-2)+iphi+1;
        GN4=GNleft+(in)*(Nphi-2)+iphi;

        Xnodes([GN1 GN2 GN3 GN4])=[Xr(in,iphi) Xr(in,iphi+1) Xr(in+1,iphi+1) Xr(in+1,iphi)];
        Ynodes([GN1 GN2 GN3 GN4])=[Yr(in,iphi) Yr(in,iphi+1) Yr(in+1,iphi+1) Yr(in+1,iphi)];
        Elem(ielem).x=Xnodes([GN1 GN2 GN3 GN4]);
        Elem(ielem).y=Ynodes([GN1 GN2 GN3 GN4]);
        Elem(ielem).GN=[GN1 GN2 GN3 GN4];
    end

end
Nelemright=ielem;
GNright=GNleft+Nl*(Nphi-2);

%center piece
for in =1:Nl-1
    
    for in2=1:Nl-1
        ielem=Nelemright+in2+(in-1)*(Nl-1);
        
        GN1=GNright+(in-1)*Nl+in2;
        GN2=GNright+(in-1)*Nl+in2+1;
        GN3=GNright+(in)*Nl+in2+1;
        GN4=GNright+(in)*Nl+in2;

        Xnodes([GN1 GN2 GN3 GN4])=[Xc(in,in2) Xc(in,in2+1) Xc(in+1,in2+1) Xc(in+1,in2)];
        Ynodes([GN1 GN2 GN3 GN4])=[Yc(in,in2) Yc(in,in2+1) Yc(in+1,in2+1) Yc(in+1,in2)];
        Elem(ielem).x=Xnodes([GN1 GN2 GN3 GN4]);
        Elem(ielem).y=Ynodes([GN1 GN2 GN3 GN4]);
        Elem(ielem).GN=[GN1 GN2 GN3 GN4];
    end

end
Nelemcenter=ielem;
GNcenter=GNright+Nl*Nl;

%boundary between the zones
%top left
for in=1:Nl-1
    ielem=Nelemcenter+in;
    
    GN1=(in-1)*(Nphi-2)+1;
    GN2=(in)*(Nphi-2)+1;
    GN3=GNright+in+1;
    GN4=GNright+in;

    Elem(ielem).x=Xnodes([GN1 GN2 GN3 GN4]);
    Elem(ielem).y=Ynodes([GN1 GN2 GN3 GN4]);
    Elem(ielem).GN=[GN1 GN2 GN3 GN4];
end

%bottom left
for in=1:Nl-1
    ielem=Nelemcenter+1*(Nl-1)+in;
    
    GN1=(in)*(Nphi-2);
    GN2=(in+1)*(Nphi-2);
    GN3=GNcenter-in*Nl;
    GN4=GNcenter-(in-1)*Nl;

    Elem(ielem).x=Xnodes([GN1 GN2 GN3 GN4]);
    Elem(ielem).y=Ynodes([GN1 GN2 GN3 GN4]);
    Elem(ielem).GN=[GN1 GN2 GN3 GN4];
end

%top right
for in=1:Nl-1
    ielem=Nelemcenter+2*(Nl-1)+in;
    
    GN1=GNright+Nl*(in-1)+1;
    GN2=GNright+Nl*in+1;
    GN3=GNright-(Nphi-2)*(Nl-in-1);
    GN4=GNright-(Nphi-2)*(Nl-in);

    Elem(ielem).x=Xnodes([GN1 GN2 GN3 GN4]);
    Elem(ielem).y=Ynodes([GN1 GN2 GN3 GN4]);
    Elem(ielem).GN=[GN1 GN2 GN3 GN4];
end

%bottom right
for in=1:Nl-1
    ielem=Nelemcenter+3*(Nl-1)+in;
    
    GN1=GNcenter-(in-1);
    GN2=GNcenter-in;
    GN3=GNleft+1+(Nphi-2)*(in);
    GN4=GNleft+1+(Nphi-2)*(in-1);

    Elem(ielem).x=Xnodes([GN1 GN2 GN3 GN4]);
    Elem(ielem).y=Ynodes([GN1 GN2 GN3 GN4]);
    Elem(ielem).GN=[GN1 GN2 GN3 GN4];
end


end

