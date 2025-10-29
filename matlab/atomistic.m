%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        PROGRAM FOR RESIDUAL INTERACTION PROBABILITY CALCULATION
%               USING ATOMISTIC VAN DER WAALS POTENTIALS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input files:
% - vdw_cns.param, parameter file used in MD (for obtaining atom types)
% - pdb_cns.top, topology file used in MD (for obtaining energy values)
% Output of the program:
% - normalized energy matrix
% - total energy matrix
% weight: e^-E/kt
% normalized: normalized version of totalenergies matrix

function y = atomistic(ProteinName)
File = 'vdw_cns.param';
rcut = 5; % cutoff radius (Å)

fid = fopen(File);
i = 0;
while feof(fid) == 0
    tline = fgets(fid,255);
    if  tline(2:10) == 'NONBonded'
        i = i+1;
        partype(i,:)  = tline(13:16);
        pareps(i,1)   = str2num(tline(21:26));
        parsigma(i,1) = str2num(tline(30:35));
    end
end
fclose(fid);

File = 'pdb_cns.top';
fid = fopen(File);
i = 0;
while feof(fid) == 0
    tline = fgets(fid,255);
    i = i+1;
    parrestype(i,:) = tline(1:3);
    paratomtype(i,:) = tline(11:13);
    partypecor(i,:) = tline(21:24); % corresponding type to atomtype
    charge(i,:) = tline(34:38);
end
fclose(fid); 

kT = 1;
pdb = pdbreader(ProteinName);
k = 1;
mult_count = 0;
for i = 1:length(pdb.Model.Atom)
    if strcmp(pdb.Model.Atom(1,i).iCode,'') & strcmp(pdb.Model.Atom(1,i).altLoc,'')
        atomnumber1(k,1) = pdb.Model.Atom(1,i).AtomSerNo;
        atomtype1{k}     = upper(pdb.Model.Atom(1,i).AtomName);
        altloc{k} = [];
        if length(pdb.Model.Atom(1,i).altLoc) > 0
            altloc{k} = pdb.Model.Atom(1,i).altLoc;
        end
        restype1{k}     = upper(pdb.Model.Atom(1,i).resName);
        chainid1{k}     = upper(pdb.Model.Atom(1,i).chainID);
        resnumber1(k,1) = pdb.Model.Atom(1,i).resSeq;
        x1(k,1)         = pdb.Model.Atom(1,i).X;
        y1(k,1)         = pdb.Model.Atom(1,i).Y;
        z1(k,1)         = pdb.Model.Atom(1,i).Z;
        bfactor1(k,1)   = 1;
        rfactor1(k,1)   = 1;
        element1{i}     = pdb.Model.Atom(1,i).element;
        k = k+1;
    else
        mult_count = mult_count+1;
    end
end
ds = length(pdb.Model.Atom) - mult_count;

CorFile = strcat(ProteinName, '.cor');
rcor = dlmread(CorFile); % input obtained from readpdb.f

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filter hydrogens/altLocs as original logic
sum  = 0; sumx = 0; sumy = 0; sumz = 0; sumrf = 0;
j = 0; k = -1; i = 0;
while i ~= ds
    i = i+1;
    j = j+1;
    if strcmp(restype1{i},'HIP')==1 || strcmp(restype1{i},'HIE')==1 || strcmp(restype1{i},'HID')==1 || strcmp(restype1{i},'HSD')==1 || strcmp(restype1{i},'HSE')==1 || strcmp(restype1{i},'HSP')==1 
        restype1{i} = 'HIS';
    elseif strcmp(restype1{i},'LYN')==1
        restype1{i} = 'LYS';
    elseif strcmp(restype1{i},'GLH')==1 
        restype1{i} = 'GLU';
    elseif strcmp(restype1{i},'CYM')==1 || strcmp(restype1{i},'CYX')==1
        restype1{i} = 'CYS';
    elseif strcmp(restype1{i},'ASH')==1
        restype1{i} = 'ASN';
    elseif strcmp(restype1{i},'TYM')==1
        restype1{i} = 'TYR';
    end
    if i ~= ds
        while strcmp(element1{i},'H')==1 || strcmp(atomtype1{i},'H')==1 || strcmp(atomtype1{i},'D')==1 || strcmp(atomtype1{i},'E')==1 || strcmp(atomtype1{i},'G')==1 || strcmp(atomtype1{i},'HN')==1 || (any(atomtype1{i}=='H')==1 & ~any(atomtype1{i}=='N'))
            i = i+1; % skip hydrogens
            if ds == i, break; end
        end
    end
    if i==ds && ( strcmp(element1{i},'H')==1 || strcmp(atomtype1{i},'H')==1 || strcmp(atomtype1{i},'D')==1 || strcmp(atomtype1{i},'E')==1 || strcmp(atomtype1{i},'G')==1 || strcmp(atomtype1{i},'HN')==1 || (any(atomtype1{i}=='H')==1 & ~any(atomtype1{i}=='N')) )
        j = j-1;
        break
    end
    if length(altloc{i})~=0 & length(altloc{i+1})~=0 & strcmp(altloc{i},altloc{i+1})==0 
        j = j-1;
        continue;
    end
    atomnumber(j,1) = atomnumber1(i,1);
    atomtype{j}     = atomtype1{i};
    restype{j}      = restype1{i};
    chainid{j}      = chainid1{i};
    resnumber(j,1)  = resnumber1(i,1);
    x(j,1)          = x1(i,1);
    y(j,1)          = y1(i,1);
    z(j,1)          = z1(i,1);
    bfactor(j,1)    = bfactor1(i,1);
    rfactor(j,1)    = rfactor1(i,1);
    element{j}      = element1{i};
end
ds = j;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Map LJ params
for i = 1:ds
    if strcmp(atomtype{i},'OXT')==1 || strcmp(atomtype{i},'O2')==1 || strcmp(atomtype{i},'O1')==1 || strcmp(element{i},'O')==1
        atomtype{i} = 'O';
    end
    for j = 1:247
        ra = safe_strtrim(parrestype(j,:));
        pa = safe_strtrim(paratomtype(j,:));
        if strcmp(restype{i}, ra)==1 & strcmp(pa, atomtype{i})==1
            for k = 1:34
                if strcmp(char(partypecor(j,:)), char(partype(k,:)))==1
                    eps(i,1)   = pareps(k,1);
                    sigma(i,1) = parsigma(k,1);
                end
            end
        end
    end
end

% Keep only atoms that received parameters (original code already resets ds)
ds = size(eps,1);

fprintf('\n\nReading epsilon and sigma values are over\n\n');
w = 0;
for i = 1:ds
    if eps(i,1)==0 || sigma(i,1)==0
        w = w+1;
        error(w,1) = atomnumber(i,1);
        berror{w}  = restype{i};
    end
end
if w > 0
    fprintf('\n\nThere is an error\n\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ENERGY CALCULATION (LJ 12-6)
for i = 1:ds
    for j = 1:ds
        r(i,j) = sqrt((x(i,1)-x(j,1))^2 + (y(i,1)-y(j,1))^2 + (z(i,1)-z(j,1))^2);
        sigmamix(i,j) = (sigma(i,1)+sigma(j,1))/2;
        epsmix(i,j)   = sqrt((eps(i,1)*eps(j,1)));
        rmin(i,j)     = (2^(1/6))*sigmamix(i,j);
        if i==j
            energy(i,j) = 0;
        elseif r(i,j) < rmin(i,j)
            energy(i,j) = 4*epsmix(i,j)*((sigmamix(i,j)/rmin(i,j))^12 - (sigmamix(i,j)/rmin(i,j))^6);
        elseif r(i,j) > rcut
            energy(i,j) = 0;
        else
            energy(i,j) = 4*epsmix(i,j)*((sigmamix(i,j)/r(i,j))^12 - (sigmamix(i,j)/r(i,j))^6);
        end
    end
    if mod(i,1000)==0
        fprintf('Atoms calculated so far: %d\n\n', i);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RESIDUAL ENERGY CALCULATION
A = energy;                 % (ds x ds)
[Asz1, Asz2] = size(A);     % should be ds x ds
xsz = size(rcor);
totresno = xsz(1);

resenergy = zeros(totresno, totresno); % preallocate
sum_val = 0;
kalani = 1;
kalanj = 1;

for l = 1:totresno
    for k = 1:totresno
        % Intended atom-block indices from .cor
        i_start = kalani;
        i_end   = kalani + rcor(l,9) - 1;
        j_start = kalanj;
        j_end   = kalanj + rcor(k,9) - 1;

        % ---- Bounds clamp to available matrix A (minimal change, avoids OOB)
        if i_start > Asz1 || j_start > Asz2
            sum_val = 0; % nothing to add if block starts beyond A
        else
            i_end_eff = min(i_end, Asz1);
            j_end_eff = min(j_end, Asz2);
            sum_val = 0;
            for ii = i_start:i_end_eff
                for jj = j_start:j_end_eff
                    sum_val = sum_val + A(ii,jj);
                end
            end
        end

        resenergy(l,k) = sum_val/2;
        if l == k
            resenergy(l,k) = 0;
        end

        if k ~= totresno
            kalanj = kalanj + rcor(k,9);
        end
    end
    kalanj = 1;
    kalani = kalani + rcor(l,9);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ENERGY NORMALIZATION
sum_row = 0;
for i = 1:totresno
    for j = 1:totresno
        weight(i,j) = exp(-resenergy(i,j)/kT);  % weight = exp(-E/kt)
        sum_row = weight(i,j) + sum_row;        % row sum
    end
    for j = 1:totresno
        normalized(i,j) = weight(i,j)/sum_row;  % normalization
        if i == j
            normalized(i,j) = 0;
        end
    end
    sum_row = 0;
end

File = strcat(ProteinName, '_atomistic.out');
dlmwrite(File, normalized, '\t'); % probability matrix
y = 0; % keep original function signature
end

% ------------------------- Helper (Octave-safe) --------------------------
function s = safe_strtrim(x)
% Accept char, string, or numeric arrays and return a trimmed char row.
  if ischar(x)
    s = strtrim(x);
  elseif isstring(x)
    s = strtrim(char(x));
  elseif isnumeric(x)
    s = strtrim(char(x));
  else
    s = strtrim(char(x));
  end
end
