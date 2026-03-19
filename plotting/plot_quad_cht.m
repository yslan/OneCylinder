function f=plot_quad(ifig,X,Quad,CBC_list)
if(ifig==0);return;end; if(ifig<0); ifig=abs(ifig); hold off; end; f=figure(ifig);

   CBCv=CBC_list{1}; nelgv = size(CBCv,1);
   CBCt=CBC_list{2}; nelgt = size(CBCt,1);

   nX=size(X,1); dim = size(X,2);

   id=min((Quad>0) & (Quad<=nX),[],2);
   Quad = Quad(id,:);

   nQ = nelgv;
   iq = 1:nelgv;
   Xq = reshape(X(Quad(iq,[1,2,3,4,1]),1),nQ,5);
   Yq = reshape(X(Quad(iq,[1,2,3,4,1]),2),nQ,5);

   fill(Xq',Yq','w'); hold on
   plot(Xq',Yq','k-'); hold on

   nQ = nelgt - nelgv;
   iq = (nelgv+1):nelgt;
   Xq = reshape(X(Quad(iq,[1,2,3,4,1]),1),nQ,5);
   Yq = reshape(X(Quad(iq,[1,2,3,4,1]),2),nQ,5);
   fill(Xq',Yq','r'); hold on
   plot(Xq',Yq','k-'); hold on

   axis equal
   grid on
