warning off;clear all; close all; format compact; profile off; diary off; restoredefaultpath;warning on;
pause(.1);hdr;

verbose = 1;
ifgui = 0;
ifpng = 1; % save pngs

% debug
fldr_out = 'outputs_cyl2d';

% fluid
   nbx0 = 10; % E_box = nbx0 x nbx0
L0 = 0.4; % inner box: [-L, L]
   nlayers1 = 6;
   ratio1 = 1/1.1;
R1 = 0.9; % bdry layers
   nlayers2 = 3;
   ratio2 = 1/1.2;
R2 = 0.99; % bdry layer
   nlayers3 = 1;
   ratio3 = 1;
Rcyl = 1.0;

% solid
   nlayers4 = 3;
   ratio4 = 1.3;
Rt = 1.1; % thickness



% for thanh1p
%Rt = 1.0 + 0.05637707948;
%fldr_out = 'outputs_cyl2d_reso2';

%   nbx0 = 6; % E_box = nbx0 x nbx0
%   nlayers1 = 3;
%   nlayers2 = 1;
%R2 = 0.97;
%Rt = 1.0 + 0.05637707948;
%fldr_out = 'outputs_cyl2d_reso1';


mkdir(fldr_out);
if(ifgui==0);set(0,'DefaultFigureVisible','off');end % server mode

save_aux = @(h,tag) print(h, '-dpng','-r400',[fldr_out '/' tag '.png']);
if (ifpng==0); save_aux = @(h,tag)[]; end

% generate box mesh
tag = 'box'; deform = 0.45 * 1/sqrt(2); % curvature of the box: 1/sqrt(2) = circle
[X, Quad, Qfront] = gen_box2d(nbx0, deform); X = X * L0; nQ = size(Quad,1); 
Qcurve = zeros(6, 4, nQ); % (type + bc(5), 4 faces, E)
Qbc = zeros(4, nQ);
bc_set = chk_bcid([],Qbc,tag,1);
ifig = 1; h=plot_quad(ifig,X,Quad); save_aux(h,tag);

% box to cyl
tag = 'cir1'; Nlap = 0;
[X, Quad, Qbc, Qcurve, Qfront] = extrude2d_front2cyl(X, Quad, Qbc, Qcurve, Qfront, R1, nlayers1, ratio1, Nlap);
bc_set = chk_bcid([],Qbc,tag,1);
ifig = 2; h=plot_quad(ifig,X,Quad); save_aux(h,tag);
nQ1 = size(Quad,1);

% circle2
tag = 'cir2'; Nlap = 0;
[X, Quad, Qbc, Qcurve, Qfront] = extrude2d_front2cyl(X, Quad, Qbc, Qcurve, Qfront, R2, nlayers2, ratio2, Nlap);
bc_set = chk_bcid([],Qbc,tag,1);
ifig = 3; h=plot_quad(ifig,X,Quad); save_aux(h,tag);

% to bdry
tag = 'cir3'; Nlap = 0; 
[X, Quad, Qbc, Qcurve, Qfront] = extrude2d_front2cyl(X, Quad, Qbc, Qcurve, Qfront, Rcyl, nlayers3, ratio3, Nlap);
bc_set = chk_bcid([],Qbc,tag,1);
ifig = 4; h=plot_quad(ifig,X,Quad); save_aux(h,tag);

% maps faces to CBC
nQ = size(Quad,1);
nelv = size(Quad,1);

bc_map = {};
bc_map{1} = [9,10,11,12];
CBCv = set_cbc(Qbc,bc_map)';

% solid domain
tag = 'cir4'; Nlap = 0; 
[X, Quad, Qbc, Qcurve, Qfront] = extrude2d_front2cyl(X, Quad, Qbc, Qcurve, Qfront, Rt, nlayers4, ratio4, Nlap);
bc_set = chk_bcid([],Qbc,tag,1);
ifig = 5; h=plot_quad(ifig,X,Quad); save_aux(h,tag);

bc_map = {};
bc_map{1} = [13,14,15,16];
CBCt = set_cbc(Qbc,bc_map)';

% reconstruct cyl curves for face 1 and 3
nQ = size(Quad,1);
iq = (nQ1+1):nQ;
Qcurve_tmp = fix_cyl_curve(X,Quad(iq,:),Qcurve(:,:,iq));
Qcurve(:,:,iq) = Qcurve_tmp;



% TODO: n2to3

%% dump mesh
fout='out';
fname=[fldr_out '/' fout];
dump_nek_rea_heat(fname,X,Quad,{CBCv,CBCt},1,Qcurve,verbose+1);

flist = {'nekwriter/head.rea', [fname '.out'], 'nekwriter/tail.rea'}; % TODO chk 3D head.rea, ifheat, nBC, etc
fout = [fname '.rea'];
combine_txt(flist, fout, verbose)


