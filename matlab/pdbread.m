% pdbread.m (Octave shim) — minimal ATOM/HETATM parser
% Returns struct pdb with pdb.Model.Atom(k).fields used by your code:
% AtomSerNo, AtomName, altLoc, resName, chainID, resSeq, X, Y, Z, element, iCode
function pdb = pdbread(filename)
  fid = fopen(filename, 'r');
  assert(fid > 0, 'pdbread: cannot open file %s', filename);
  Atoms = struct('AtomSerNo',{},'AtomName',{},'altLoc',{},'resName',{}, ...
                 'chainID',{},'resSeq',{},'X',{},'Y',{},'Z',{}, ...
                 'element',{},'iCode',{});
  k = 0;
  while true
    t = fgetl(fid);
    if ~ischar(t), break; end
    if length(t) >= 54 && (strncmp(t,'ATOM  ',6) || strncmp(t,'HETATM',6))
      k = k + 1;
      % Safe substring helpers
      s = @(a,b) strtrim(t(a:min(b, length(t))));
      Atoms(k).AtomSerNo = str2double(s(7,11));
      Atoms(k).AtomName  = s(13,16);
      Atoms(k).altLoc    = s(17,17);
      Atoms(k).resName   = s(18,20);
      Atoms(k).chainID   = s(22,22);
      Atoms(k).resSeq    = str2double(s(23,26));
      Atoms(k).iCode     = s(27,27);
      Atoms(k).X         = str2double(s(31,38));
      Atoms(k).Y         = str2double(s(39,46));
      Atoms(k).Z         = str2double(s(47,54));
      % element: try columns 77–78 else derive from atom name
      if length(t) >= 78
        elem = s(77,78);
      else
        nm = regexprep(Atoms(k).AtomName,'[^A-Za-z]','');
        elem = ''; if ~isempty(nm), elem = nm(1); end
      end
      Atoms(k).element = elem;
    end
  end
  fclose(fid);
  pdb.Header = struct();           % not used by your script
  pdb.Model  = struct('Atom', Atoms);
end
