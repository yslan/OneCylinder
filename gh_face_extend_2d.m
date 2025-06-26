function x = gh_face_extend_2d(x,zg,n,gh_type)
% copy from the subroutine gh_face_extend_2d from Nek5000 / navier5.f 
%
%     Extend 2D faces into interior via gordon hall
%
%     gh_type:  1 - vertex only
%               2 - vertex and faces

   sz = size(x);
   x = reshape(x,n,n);
   v = zeros(n,n);      
%   [zg,~] = zwgll(n-1);

%
%     Build vertex interpolant
%
   for jj=[1,n];
   for ii=[1,n];
      for j=1:n
      for i=1:n
            si     = 0.5*((n-ii)*(1-zg(i))+(ii-1)*(1+zg(i)))/(n-1);
            sj     = 0.5*((n-jj)*(1-zg(j))+(jj-1)*(1+zg(j)))/(n-1);
            v(i,j) = v(i,j) + si*sj*x(ii,jj);
      end
      end
   end
   end
   
   if (gh_type==1);
      x = v;
      x = reshape(x,sz);
      return
   end


%     Extend 4 edges
%     call rzero(e,ntot)
   e = zeros(n,n);
%
%     x-edges
%     
   for jj=[1,n]
      for j=1:n
      for i=1:n
            hj     = 0.5*((n-jj)*(1-zg(j))+(jj-1)*(1+zg(j)))/(n-1);
            e(i,j) = e(i,j) + hj*(x(i,jj)-v(i,jj));
      end
      end
   end
%
%     y-edges
%     
   for ii=[1,n]
      for j=1:n
      for i=1:n
            hi     = 0.5*((n-ii)*(1-zg(i))+(ii-1)*(1+zg(i)))/(n-1);
            e(i,j) = e(i,j) + hi*(x(ii,j)-v(ii,j));
      end
      end
   end
   x = e + v;      
   x = reshape(x,sz);


