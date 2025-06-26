function [X, Quad, Qbc, Qcurve, Qfront] = extrude2d_front2cyl(X, Quad, Qbc, Qcurve, Qfront, Rcyl, nlayers, ratio, Nlap);

if (nlayers==0); return; end

fpts_front = Qfront.fpts_front;
eids = Qfront.eids;
fids = Qfront.fids;

nptf = size(fpts_front,1);
nelf = nptf - 1;
nface = size(fpts_front,2);

nQ0 = size(Quad,1);
nQnew = nelf * nface * nlayers;

fprintf('ext2cyl: nQ= %d -> %d, nlyr= %d', nQ0, nQ0 + nQnew, nlayers);
if (length(ratio)==1)
   wt = zeros(nlayers,1);
   wt(1) = 1.0;
   for i=2:nlayers
      wt(i) = ratio * wt(i-1);
   end
   wt = wt / sum(wt);
   
   fprintf(', interp: %2.2f',wt(1));
   wt_acc = zeros(nlayers,1);
   wt_acc(1) = wt(1);
   for i=2:nlayers
      wt_acc(i) = wt_acc(i-1) + wt(i);
      fprintf(' %2.2f',wt_acc(i));
   end
   fprintf(' \n');
else % custom ratio
   wt_acc = ratio;
   assert(nlayers==length(wt_acc), 'extrude2d_box2cyl: invalid custom ratio');
   fprintf(', interp c: %2.2f',wt_acc(1));
   assert(wt_acc(1) > 0, 'extrude2d_box2cyl: invalid custom ratio');
   for i=2:nlayers
      assert(wt_acc(i) > wt_acc(i-1), 'extrude2d_box2cyl: invalid custom ratio');
      fprintf(' %2.2f',wt_acc(i));
   end
   fprintf(' \n');
end
wt_acc(end) = 1.0; 

Qnew = zeros(nQnew,4);iq0=0;
Qcurve_new = zeros(6, 4, nQnew); % (type + bc(5), 4 faces, E)
Qbc_new = zeros(4, nQnew);

tmp = abs(Qbc(:)); tmp=tmp(tmp>0); nbc0 = length(unique(tmp));

nX0 = size(X,1);
nXnew = nptf * nface * nlayers;
Xnew = zeros(nXnew,2);ix0=0;

fpts_front_new = zeros(nptf, nface);
eids_new = zeros(nelf, nface);
fids_new = zeros(nelf, nface);

for f = 1:nface
   vf = fpts_front(:,f);
   Xfbot = X(vf,:);
%   Rad = sqrt(sum(Xfbot.^2,2));
%   Xftop = Xfbot * Rcyl ./ Rad;
   Xftop = gen_x_on_cyl(Rcyl, Xfbot(1,:), Xfbot(end,:), nelf);

   v1=vf(1:nelf);
   v2=vf(2:nelf+1);
   for ilyr = 1:nlayers
      w = wt_acc(ilyr);
      Xnew(ix0+1:ix0+nptf,:) = Xfbot * (1.0 - w) + Xftop * w;

      vx = nX0 + (ix0+1:ix0+nptf)';
      v3 = vx(2:end);
      v4 = vx(1:end-1);
      ix0 = ix0 + nptf;
      
      Qnew(iq0+1:iq0+nelf,:) = [v1, v2, v3, v4];
      if (ilyr == 1); % copy curve, TODO: refactor this
         for ef = 1:nelf
           ee = eids(ef,f);
           ff = fids(ef,f);
           ctype = Qcurve(1,ff,ee);
           if Qcurve(1,ff,ee) > 0
             Qcurve_new(:,1,iq0+ef) = Qcurve(:,ff,ee);
           end
           if ctype==2; % cylinder, concave
             Qcurve_new(2,1,iq0+ef) = -abs(Qcurve(2,ff,ee));
           end
         end
      end
      if (ilyr == nlayers);
         Qbc_new(3, iq0+1:iq0+nelf) = nbc0 + f;
         Qcurve_new(1, 3, iq0+1:iq0+nelf) = 2;
         Qcurve_new(2, 3, iq0+1:iq0+nelf) = Rcyl;

         eids_new(:,f) = nQ0 + (iq0+1:iq0+nelf);
         fids_new(:,f) = 3;
      end
      v1 = v4;
      v2 = v3;
      iq0 = iq0 + nelf;
   end

   fpts_front_new(:,f) = vx;
end

assert(ix0 == nXnew, 'bad estimation of nXnew');
assert(nXnew == size(Xnew,1), 'bad estimation of nXnew');
assert(iq0 == nQnew, 'bad estimation of nQnew');
assert(nQnew == size(Qnew,1), 'bad estimation of nQnew');

X = [X;Xnew];
Quad = [Quad; Qnew];
Qbc = cat(2, Qbc, Qbc_new);
Qcurve = cat(3, Qcurve, Qcurve_new); 

Qfront.fpts_front = fpts_front_new;
Qfront.eids = eids_new;
Qfront.fids = fids_new;


X = smooth_block_quad(X, Qnew, Nlap);


function X = gen_x_on_cyl(Rcyl, x1, x2, nelf)

nptf = nelf + 1;

% equi-angle
theta = pi/2 / nelf;
r = @(i) sin(i* theta) ./ sin(3*pi/4-i*theta);
wt_theta = r(0:nelf) / sqrt(2);

xd = x2 - x1;
X = zeros(nptf, 2);
for i=1:nptf
   X(i,:) = x1 + xd*wt_theta(i);
end

R = sqrt(sum(X.^2,2));
X = X * Rcyl ./ R;


function X = smooth_block_quad(X0, quad, Nlap)

persistent ifcalld J PERM mask wght nq0

nq=size(quad,1); if(Nlap==0);X=X0;return;end

if (isempty(ifcalld) || ifcalld==0 || nq~=nq0); ifcalld=1; nq0=nq;
  perm=[1 2 4 3]; PERM=eye(4); PERM=PERM(:,perm);

% xm=-0.95; xp=0.95;
  xm=-0.9; xp=0.9;
  Jh=zeros(2,2);

  Jh(1,1)=.5*(1-xm); Jh(1,2)=.5*(1+xm);
  Jh(2,1)=.5*(1-xp); Jh(2,2)=.5*(1+xp);

  J=PERM* kron(Jh,Jh) *PERM;

%  mask=1+0*quad; mask=gs_qqt_mask(mask,quad,edge,1);
  mask=1+0*quad; mask=gs_qqt(mask,quad);
  jd0=mask==0; wght=1./mask;wght(jd0)=0;
  jd1=mask<=2; jd2=mask>2; mask(jd1<=2)=0;mask(jd2)=1; % mask: bdry=0, int=1
end

nQ=size(quad,1); dim=size(X0,2);
X=X0;for d=1:dim;Xtmp=X0(quad,d);Xl0(:,:,d)=reshape(Xtmp,nQ,4);end

for d=1:dim; iprog=1;fprintf('    adjQ, smooth dim-%d, progressing ...',d);tp=tic;

  Xtmp=X(quad,d); Xl(:,:,d)=reshape(Xtmp,nQ,4); % GtoL
  for ilap=1:Nlap
    Xd=Xl(:,:,d);
    Xd=Xd*J';
    Xd=wght.*gs_qqt(Xd,quad);
    Xl(:,:,d)=Xd.*mask + (1-mask).*Xl0(:,:,d);
    if(ilap/Nlap>=iprog/10);fprintf(' %d',iprog);iprog=iprog+1;end
  end; fprintf(' done! %2.4e sec\n',toc(tp));

  X(quad,d)=reshape(Xl(:,:,d),nQ*4,1); % LtoG
end
