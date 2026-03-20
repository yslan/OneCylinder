warning off;clear all; close all; format compact; profile off; diary off; restoredefaultpath;warning on;
pause(.1);hdr;

verbose = 1;
ifgui = 1;
ifpng = 1; % save pngs

% debug
fldr_out = 'outputs_double_cyl2d';

% fluid
   nbx0 = 4; % E_box = nbx0 x nbx0
L0 = 0.4; % inner box: [-L, L]
   nlayers1 = 2;
   ratio1 = 1/1.3;
R1 = 0.9; % bdry layers
   nlayers2 = 2;
   ratio2 = 1/1.3;
R2 = 1.0;

% solid
   nlayers3 = 3;
   ratio3 = 1.0;
R3 = 1.2; % thickness

% second fluid
   nlayers4 = 4;
   ratio4 = 1.6;
R4 = 1.6;
   nlayers5 = 4;
   ratio5 = 1/1.6;
R5 = 2.0;

mkdir(fldr_out);
if(ifgui==0);set(0,'DefaultFigureVisible','off');end % server mode

apply_xylim = @(x) axis([-2.2 2.2 -2.2 2.2]);

save_aux = @(h,tag) print(h, '-dpng','-r400',[fldr_out '/' tag '.png']);
if (ifpng==0); save_aux = @(h,tag)[]; end

% generate box mesh
tag = 'box'; deform = 0.45 * 1/sqrt(2); % curvature of the box: 1/sqrt(2) = circle
[X, Quad, Qfront] = gen_box2d(nbx0, deform); X = X * L0; nQ = size(Quad,1); 
Qcurve = zeros(6, 4, nQ); % (type + bc(5), 4 faces, E)
Qbc = zeros(4, nQ);
bc_set = chk_bcid([],Qbc,tag,1);
ifig = 1; h=plot_quad(ifig,X,Quad); apply_xylim(); save_aux(h,tag);

% box to cyl
tag = 'cir1'; Nlap = 0;
[X, Quad, Qbc, Qcurve, Qfront] = extrude2d_front2cyl(X, Quad, Qbc, Qcurve, Qfront, R1, nlayers1, ratio1, Nlap);
bc_set = chk_bcid([],Qbc,tag,1);
ifig = 2; h=plot_quad(ifig,X,Quad); apply_xylim(); save_aux(h,tag);
nQ1 = size(Quad,1); % for fixing curves

% to bdry
tag = 'cir2'; Nlap = 0; 
[X, Quad, Qbc, Qcurve, Qfront] = extrude2d_front2cyl(X, Quad, Qbc, Qcurve, Qfront, R2, nlayers2, ratio2, Nlap);
bc_set = chk_bcid([],Qbc,tag,1);
ifig = 3; h=plot_quad(ifig,X,Quad); apply_xylim(); save_aux(h,tag);

nel1 = size(Quad,1);


% solid domain
tag = 'cir3'; Nlap = 0; 
[X, Quad, Qbc, Qcurve, Qfront] = extrude2d_front2cyl(X, Quad, Qbc, Qcurve, Qfront, R3, nlayers3, ratio3, Nlap);
bc_set = chk_bcid([],Qbc,tag,1);
ifig = 4; h=plot_quad(ifig,X,Quad); apply_xylim(); save_aux(h,tag);

nel2 = size(Quad,1); 

% second fluid
tag = 'cir4'; Nlap = 0; 
[X, Quad, Qbc, Qcurve, Qfront] = extrude2d_front2cyl(X, Quad, Qbc, Qcurve, Qfront, R4, nlayers4, ratio4, Nlap);
bc_set = chk_bcid([],Qbc,tag,1);
ifig = 5; h=plot_quad(ifig,X,Quad); apply_xylim(); save_aux(h,tag);

tag = 'cir5'; Nlap = 0; 
[X, Quad, Qbc, Qcurve, Qfront] = extrude2d_front2cyl(X, Quad, Qbc, Qcurve, Qfront, R5, nlayers5, ratio5, Nlap);
bc_set = chk_bcid([],Qbc,tag,1);
ifig = 6; h=plot_quad(ifig,X,Quad); apply_xylim(); save_aux(h,tag);

nel3 = size(Quad,1); 

bc_map = {};
bc_map{1} = [5,6,7,8];
bc_map{2} = [-13,-14,-15,-16];
bc_map{3} = [17,18,19,20];
CBCv = set_cbc(Qbc,bc_map)';

bc_map = {};
bc_map{1} = [17,18,19,20];
CBCt = set_cbc(Qbc,bc_map)';

% reconstruct cyl curves for face 1 and 3
nQ = size(Quad,1);
iq = (nQ1+1):nQ;
Qcurve_tmp = fix_cyl_curve(X,Quad(iq,:),Qcurve(:,:,iq));
Qcurve(:,:,iq) = Qcurve_tmp;


% reorder elements so solid is placed last
iq = [(1:nel1), (nel2+1):nel3, (nel1+1):nel2];
nels = nel2 - nel1;
nelt = nel3;
nelv = nelt - nels;
CBCv = CBCv(iq,:); CBCv = CBCv(1:nelv,:);
CBCt = CBCt(iq,:);

Qcurve = Qcurve(:,:,iq);
Quad = Quad(iq,:);
ifig = 50; tag='domain';
   h=plot_quad_cht(ifig,X,Quad,{CBCv,CBCt}); apply_xylim(); save_aux(h,tag);

ifig = 60; tag='fluid_bc';
   h=plot_quad(ifig,X,Quad);
   h=plot_quad_bc(ifig,X,Quad,CBCv,1,'r--'); 
   h=plot_quad_bc(ifig,X,Quad,CBCv,2,'b-'); 
   h=plot_quad_bc(ifig,X,Quad,CBCv,3,'g:'); 
   apply_xylim(); save_aux(h,tag);
ifig = 61; tag='heat_bc';
   h=plot_quad(ifig,X,Quad);
   h=plot_quad_bc(ifig,X,Quad,CBCt,1,'r--'); 
   h=plot_quad_bc(ifig,X,Quad,CBCt,2,'b-'); 
   h=plot_quad_bc(ifig,X,Quad,CBCt,3,'g:'); 
   apply_xylim(); save_aux(h,tag);

% TODO: n2to3

%% dump mesh
fout='out';
fname=[fldr_out '/' fout];
dump_nek_rea_heat(fname,X,Quad,{CBCv,CBCt},1,Qcurve,verbose+1);

flist = {'nekwriter/head.rea', [fname '.out'], 'nekwriter/tail.rea'}; % TODO chk 3D head.rea, ifheat, nBC, etc
fout = [fname '.rea'];
combine_txt(flist, fout, verbose)


