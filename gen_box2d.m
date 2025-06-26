function [X,Quad,Qfront] = gen_box2d(nbx, deform)
% generate 2d uniform box mesh for [-1,1]^2, E = nbx * nbx

v  = sqrt(2)/2;

nel = nbx * nbx;
nelf = nbx;
nface = 4;

npt = nbx+1;
npts = (nbx+1)*(nbx+1);


xx = linspace(-1,1,nbx+1);
[x, y] = ndgrid(xx,xx);


X = [x(:), y(:)];


if (deform>0)
  X = deform_box(X, nelf, deform);
end

vid = reshape(1:npts,npt,npt);

v1 = vid(1:nbx,  1:nbx);
v2 = vid(2:nbx+1,1:nbx);
v3 = vid(2:nbx+1,2:nbx+1);
v4 = vid(1:nbx,  2:nbx+1);

Quad = [v1(:), v2(:), v3(:), v4(:)];


fpts_front = zeros(npt, nface); % ordered from left to right, facing outward
fpts_front(:,1) = vid(end:-1:1,1);
fpts_front(:,2) = vid(end,end:-1:1)';
fpts_front(:,3) = vid(:,end);
fpts_front(:,4) = vid(1,:)';

eids = zeros(nelf, nface);
fids = zeros(nelf, nface);
eids(:,1) = nelf:-1:1;
eids(:,2) = nelf*nelf:-nelf:nelf;
eids(:,3) = (nelf-1)*nelf+1:nelf*nelf;
eids(:,4) = 1:nelf:(nelf-1)*nelf+1;;

fids(:,1) = 1;
fids(:,2) = 2;
fids(:,3) = 3;
fids(:,4) = 4;

Qfront.fpts_front = fpts_front;
Qfront.eids = eids;
Qfront.fids = fids;


function Xnew = deform_box(X0,nelf,deform)

E = 1; nh = nelf+1;

X = X0(:,1);
Y = X0(:,2);

deform_ep = deform;
deform_wave = 0.5;
fun_deform = @(t) deform_ep * sin( deform_wave*pi * (t+1) );

%Rc = sqrt(2) * (1 + 1/deform);
Rc = 1/deform;
yc = -sqrt(Rc^2-1);
fun_deform = @(t) yc + sqrt(Rc^2-t.^2);

fprintf('deform: curvature rad= %g', Rc);

dy = fun_deform(X) .* sign(Y);
dx = fun_deform(Y) .* sign(X);

X1 = X + dx;
Y1 = Y + dy;

X1 = reshape(X1,nh,nh,E);
Y1 = reshape(Y1,nh,nh,E);

for kpass=1:3

   X = reshape(X,nh,nh,E);
   Y = reshape(Y,nh,nh,E);

   for e=1:E
      X(1,:,e)    = X1(1,:,e);
      X(nh,:,e)   = X1(nh,:,e);
      Y(:,1,e)    = Y1(:,1,e);
      Y(:,nh,e)   = Y1(:,nh,e);
   end


   X = reshape(X,nh*nh,E);
   Y = reshape(Y,nh*nh,E);

   %[zg,~]=zwgll(N);
   zg = linspace(-1,1,nh);

   for e=1:E
      X(:,e) = gh_face_extend_2d(X(:,e),zg,nh,kpass);
      Y(:,e) = gh_face_extend_2d(Y(:,e),zg,nh,kpass);
   end

end

Xnew = [X(:),Y(:)];

del = max(abs(Xnew-X0));
fprintf('  diff= (%2.4e %2.4e)\n', del);
