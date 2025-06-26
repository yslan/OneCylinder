function combine_txt(flist, fout, verbose)

if (verbose>0); fprintf('Combine files...'); end
t0=tic;

fid_out = fopen(fout, 'w');

for f = 1:length(flist);
   fname = flist{f};
   if (verbose>0); fprintf('\n  %s',fname); end
   
   fid_in = fopen(fname, 'r');
   
   if fid_in ~= -1
       data = fread(fid_in, '*char');  % read entire file as character array
       fwrite(fid_out, data);          % write directly
       fclose(fid_in);
   else
       warning('Cannot open %s', fname);
   end

end

[osize,otype]=comp_fsize(fname);
if(verbose>0); fprintf(['\ndone! (%3.1f %s %2.4e sec)\n'],osize,otype,toc(t0)); end
