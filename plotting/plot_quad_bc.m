function f=plot_quad(ifig,X,Quad,CBC,bcid,varargin)
if(ifig==0);return;end; if(ifig<0); ifig=abs(ifig); hold off; end; f=figure(ifig);

lt = 'r--'; if (length(varargin)==1); lt = varargin{1}; end

   id = find(CBC==bcid);
   nedge = length(id);
   [idq,idf] = ind2sub(size(CBC), id);

   iftoiv = [1,2;2,3;3,4;4,1];
   idv1 = iftoiv(idf,1);
   idv2 = iftoiv(idf,2);

   nQ = size(Quad,1);
   id1 = sub2ind([nQ,4], idq, idv1);
   id2 = sub2ind([nQ,4], idq, idv2);

   Xq = reshape(X([Quad(id1),Quad(id2)],1),nedge,2);
   Yq = reshape(X([Quad(id1),Quad(id2)],2),nedge,2);

   % TODO: fancy shrinkage to plot BC inside elements rather than on edges

   plot(Xq',Yq',lt,'LineWidth',2); hold on

   s=sprintf('bcid= %d  #bc= %d',bcid,nedge);
   title(s)
   axis equal
   grid on
