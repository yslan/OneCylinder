function Qcurve = fix_cyl_curve(X,Quad,Qcurve);

% face 1
tol = 1e-8;

v1 = Quad(:,1);
v2 = Quad(:,2);
v3 = Quad(:,3);
v4 = Quad(:,4);

x1 = X(v1,:);
x2 = X(v2,:);
x3 = X(v3,:);
x4 = X(v4,:);

dim = size(X,2);
E = size(Quad,1);
nv = size(Quad,2);
Xloc = zeros(E,nv,dim);;

for d=1:dim;
   Xtmp=X(Quad,d);
   Xloc(:,:,d)=reshape(Xtmp,E,nv);
end

Rad = sqrt(Xloc(:,:,1).^2 + Xloc(:,:,2).^2);
nface1 = 0;
nface3 = 0;
for e=1:E
   % face 1
   if abs(Rad(e,1) - Rad(e,2)) < tol;
      r = (Rad(e,1) + Rad(e,2)) / 2;
      Qcurve(1,1,e) = 2;
      Qcurve(2,1,e) = -r;
      nface1 = nface1 + 1;
   end

   % face 3
   if abs(Rad(e,3) - Rad(e,4)) < tol;
      r = (Rad(e,3) + Rad(e,4)) / 2;
      Qcurve(1,3,e) = 2;
      Qcurve(2,3,e) = r;
      nface3 = nface3 + 1;
   end
end

fprintf('fix_cyl_curve f1=%d f3=%d\n',nface1,nface3);
