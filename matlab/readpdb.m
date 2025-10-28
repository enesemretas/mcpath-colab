%%%%%%%%%%%%%%% READ PDB %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For Extracting Data for Calculations from PDB Files
% Input files: pdb code
% Output files:
% - Cor file - Coordinates averaged for alpha carbons
% - Cbcor file - Coordinates averaged for beta carbons
% - atomcor file - Coordinates of all atoms
% - fi file - Rotational angles

function y=readpdb(ProteinName,chainID)
y=0;
if nargin<2, chainID = []; end
%fprintf('Program for Extracting Data from PDB Files\n');
%ProteinName=input('Enter the name of the pdb file as "xxxx":\n','s');
%File = ['../upload/' ProteinName];
pdb = pdbreader(ProteinName);
    % Collect chain IDs from ATOM records (Octave-safe)
atoms = pdb.Model.Atom;
atomChains = arrayfun(@(a) upper(a.chainID), atoms, 'UniformOutput', false);  % cellstr
chains = unique(atomChains, 'stable');                 % keep first-seen order
chains = chains(~cellfun(@isempty, strtrim(chains)));  % drop blanks
rank = find(strcmpi(chains, chainID), 1);

k=1;
mult_count=0;
for i=1:length(pdb.Model.Atom)
    if strcmp(pdb.Model.Atom(1,i).iCode,'') & strcmp(pdb.Model.Atom(1,i).altLoc,'')
    atomnumber1(k,1)=pdb.Model.Atom(1,i).AtomSerNo;
    atomtype1{k}=upper(pdb.Model.Atom(1,i).AtomName);
    altloc{k}=[];
    if length(pdb.Model.Atom(1,i).altLoc)>0
        altloc{k}=pdb.Model.Atom(1,i).altLoc;
    end
    restype1{k}=upper(pdb.Model.Atom(1,i).resName);
    chainid1{k}=upper(pdb.Model.Atom(1,i).chainID);
    resnumber1(k,1)=pdb.Model.Atom(1,i).resSeq;
    x1(k,1)=pdb.Model.Atom(1,i).X;
    y1(k,1)=pdb.Model.Atom(1,i).Y;
    z1(k,1)=pdb.Model.Atom(1,i).Z;
    bfactor1(k,1)=1;
    rfactor1(k,1)=1;
    element1{k}=pdb.Model.Atom(1,i).element;
    k=k+1;
    else
        mult_count=mult_count+1;
end
ds=length(pdb.Model.Atom)-mult_count;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NEW: For averaging the coordinates for the atoms with two coordinates
sum=0;
sumx=0;
sumy=0;
sumz=0;
sumrf=0;
j=0;
k=-1;
i=0;
while i~=ds
    i=i+1;
    j=j+1;
    if strcmp(restype1{i},'HIP')==1 || strcmp(restype1{i},'HIE')==1 || strcmp(restype1{i},'HID')==1 || strcmp(restype1{i},'HSD')==1 || strcmp(restype1{i},'HSE')==1 || strcmp(restype1{i},'HSP')==1
        restype1{i}='HIS';
    elseif strcmp(restype1{i},'LYN')==1
        restype1{i}='LYS';
    elseif strcmp(restype1{i},'GLH')==1 
        restype1{i}='GLU';
    elseif strcmp(restype1{i},'CYM')==1 || strcmp(restype1{i},'CYX')==1
        restype1{i}='CYS';
    elseif strcmp(restype1{i},'ASH')==1
        restype1{i}='ASN';
    elseif strcmp(restype1{i},'TYM')==1
        restype1{i}='TYR';
    end
    if i~=ds
        while strcmp(element1{i},'H')==1 || strcmp(atomtype1{i},'H')==1 || strcmp(atomtype1{i},'D')==1 || strcmp(atomtype1{i},'E')==1 || strcmp(atomtype1{i},'G')==1 || strcmp(atomtype1{i},'HN')==1 || (any(atomtype1{i}=='H')==1 & ~any(atomtype1{i}=='N'))
            i=i+1; %for pdb codes including hydrogen atoms
            if ds==i
                break
            end
        end
    end
    if i==ds && strcmp(element1{i},'H')==1 || strcmp(atomtype1{i},'H')==1 || strcmp(atomtype1{i},'D')==1 || strcmp(atomtype1{i},'E')==1 || strcmp(atomtype1{i},'G')==1 || strcmp(atomtype1{i},'HN')==1 || (any(atomtype1{i}=='H')==1 & ~any(atomtype1{i}=='N'))
        j=j-1;
        break
    end
    if length(altloc{i})~=0 & length(altloc{i+1})~=0 & strcmp(altloc{i},altloc{i+1})==0 
        j=j-1;
        continue;
    end
%     if i~=ds & strcmp(altloc{i},' ')==0
%         if strcmp(atomtype1{i},atomtype1{i+1})==1
%             k=0;
%             while strcmp(atomtype1{i+k},atomtype1{i+k+1})==1 % if an error occurs around here the pdb file should be checked
%                 k=k+1;   % to see if the b-factors add up to be 1.00 for the same atom
%                 newb(k,1)=bfactor1(i+k-1,1);
%                 newb(k+1,1)=bfactor1(i+k,1);
%                 newx(k,1)=x1(i+k-1,1);
%                 newy(k,1)=y1(i+k-1,1);
%                 newz(k,1)=z1(i+k-1,1);
%                 newrf(k,1)=rfactor1(i+k-1,1);
%                 newx(k+1,1)=x1(i+k,1);
%                 newy(k+1,1)=y1(i+k,1);
%                 newz(k+1,1)=z1(i+k,1);
%                 newrf(k+1,1)=rfactor1(i+k,1);
%             end
%             tirt=find(newb(:,:)==max(newb));
%             tirt=tirt(1);
%             i=i+k;
%             x1(i,1)=newx(tirt,1);
%             y1(i,1)=newy(tirt,1);
%             z1(i,1)=newz(tirt,1);
%             bfactor1(i,1)=1;
%             rfactor1(i,1)=newrf(tirt,1);
%             clear newx newy newz newb newrf
%         else %%%burasi hala cok ilkel, ayni sayida atomun tamamen ayni sirayla yazilmasi lazim ki dogru sonucu versin
%             fprintf('Ben mal bir pdb yim lutfen %d nolu satirlardan itibaren kontrol edin\n\n',i);
%             k=0; %new counter
%             p=1; %new counter
%             w=0; %new counter
%             t=i; %to remember i
%             while i~=ds &&strcmp(altloc{i+p-1},' ')==0 && resnumber1(i+p,1)==resnumber1(i+p-1)
%                 p=p+1;
%                 if strcmp(altloc{i+p-1},altloc{i+p})==0 && resnumber1(i+p,1)==resnumber1(i+p-1)
%                     k=k+1;
%                     w(k,1)=p;
%                     newb(k,1)=bfactor1(i+p-1,1);
%                 end
%             end
%             k=k+1;
%             w(k,1)=p;
%             newb(k,1)=bfactor1(i+p-1,1);
%             k=find(newb(:,:)==max(newb));
%             k=k(1);
%             if size(w,1)==1
%                 i=i+w(end)-1;
%             else
%                 i=i+w(end-1)-1; % son grubun basina dondurmek icin
%             end
%             if k>1
%                newconst=w(k,1)-w(k-1,1);
%             else
%                 newconst=0;
%             end
%             for of=1:w(1,1) %buraya ayar cekmek lazım
%                 i=i+1;
%                 altloc{i}=' ';
%                 x1(i,1)=x1(t-1+of+newconst);
%                 y1(i,1)=y1(t-1+of+newconst);
%                 z1(i,1)=z1(t-1+of+newconst);
%                 bfactor1(i,1)=1;
%                 rfactor1(i,1)=rfactor1(t-1+of+w(k,1));
%             end
%             clear newb
%             if size(w,1)==1
%                 i=t;
%             else
%                 i=t+w(end-1); % son grubun basina dondurmek icin
%             end
%         end
%     end
    atomnumber(j,1)=atomnumber1(i,1);
    atomtype{j}=atomtype1{i};
    restype{j}=restype1{i};
    chainid{j}=chainid1{i};
    resnumber(j,1)=resnumber1(i,1);
    x(j,1)=x1(i,1);
    y(j,1)=y1(i,1);
    z(j,1)=z1(i,1);
    bfactor(j,1)=bfactor1(i,1);
    rfactor(j,1)=rfactor1(i,1);
    element{j}=element1{i};
   % atomtype2(j,1)=atomtype21(i,1);
end
ds=j;
% NEW: initialize so unknown/non-AA residues are marked 0
newrestype   = zeros(ds,1);   % 1..20 for amino acids, 0 for others (HOH, ligands, etc.)
resatomcount = zeros(ds,1);   % preallocate (you fill it later as you already do)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
chainno=0;
for i=2:ds
    a=strcmp(chainid{i-1},chainid{i});
    if a==0 % 0 means they are not the same, so the chain has changed
        chainno=chainno+1;
        chainend(chainno,1)=i;
    end
end
chainno=chainno+1;
chainend(chainno,1)=ds;
%%
gly='GLY';
ala='ALA';
val='VAL';
ile='ILE';
leu='LEU';
ser='SER';
thr='THR';
asp='ASP';
asn='ASN';
glu='GLU';
gln='GLN';
lys='LYS';
arg='ARG';
cys='CYS';
met='MET';
phe='PHE';
tyr='TYR';
trp='TRP';
his='HIS';
pro='PRO';

atomcount=0;
rescount=0;
glycount=1;
alacount=1;
valcount=1;
ilecount=1;
leucount=1;
sercount=1;
thrcount=1;
aspcount=1;
asncount=1;
glucount=1;
glncount=1;
lyscount=1;
argcount=1;
cyscount=1;
metcount=1;
phecount=1;
tyrcount=1;
trpcount=1;
hiscount=1;
procount=1;
res_num=[];
n=1;
m=1;

for i=1:ds
    newchainid(i,1)=m;
    for j=1:chainno
        if i==chainend(j,1)
            m=m+1;
        end
    end
    if strcmp(gly,restype{i})==1
        newrestype(i,1)=1;
        atomcount=atomcount+1;
        glyinfo(glycount,atomcount,1)=x(i,1);
        glyinfo(glycount,atomcount,2)=y(i,1);
        glyinfo(glycount,atomcount,3)=z(i,1);
        glyinfo(glycount,atomcount,4)=75;
        if (i==ds || resnumber(i,1)~=resnumber(i+1,1))==1
            glycount=glycount+1;
            rescount=rescount+1;
            res_num(i)=resnumber(i,1);
            resatomcount(rescount,1)=atomcount;
            atomcount=0;
        end
    end
    if strcmp(ala,restype{i})==1
        newrestype(i,1)=2;
        atomcount=atomcount+1;
        alainfo(alacount,atomcount,1)=x(i,1);
        alainfo(alacount,atomcount,2)=y(i,1);
        alainfo(alacount,atomcount,3)=z(i,1);
        alainfo(alacount,atomcount,4)=89.1;
        if (i==ds || resnumber(i,1)~=resnumber(i+1,1))==1
            alacount=alacount+1;
            rescount=rescount+1;
            resatomcount(rescount,1)=atomcount;
            res_num(i)=resnumber(i,1);
            atomcount=0;
        end
    end
    if strcmp(val,restype{i})==1
        newrestype(i,1)=3;
        atomcount=atomcount+1;
        valinfo(valcount,atomcount,1)=x(i,1);
        valinfo(valcount,atomcount,2)=y(i,1);
        valinfo(valcount,atomcount,3)=z(i,1);
        valinfo(valcount,atomcount,4)=117.2;
        if (i==ds || resnumber(i,1)~=resnumber(i+1,1))==1

            valcount=valcount+1;
            rescount=rescount+1;
            res_num(i)=resnumber(i,1);
            resatomcount(rescount,1)=atomcount;
            atomcount=0;
        end
    end
    if strcmp(ile,restype{i})==1
        newrestype(i,1)=4;
        atomcount=atomcount+1;
        ileinfo(ilecount,atomcount,1)=x(i,1);
        ileinfo(ilecount,atomcount,2)=y(i,1);
        ileinfo(ilecount,atomcount,3)=z(i,1);
        ileinfo(ilecount,atomcount,4)=131.2;
        if (i==ds || resnumber(i,1)~=resnumber(i+1,1))==1
            ilecount=ilecount+1;
            rescount=rescount+1;
            res_num(i)=resnumber(i,1);
            resatomcount(rescount,1)=atomcount;
            atomcount=0;
        end
    end
    if strcmp(leu,restype{i})==1
        newrestype(i,1)=5;
        atomcount=atomcount+1;
        leuinfo(leucount,atomcount,1)=x(i,1);
        leuinfo(leucount,atomcount,2)=y(i,1);
        leuinfo(leucount,atomcount,3)=z(i,1);
        leuinfo(leucount,atomcount,4)=131.2;
        if (i==ds || resnumber(i,1)~=resnumber(i+1,1))==1
            leucount=leucount+1;
            rescount=rescount+1;
            res_num(i)=resnumber(i,1);
            resatomcount(rescount,1)=atomcount;
            atomcount=0;
        end
    end
    if strcmp(ser,restype{i})==1
        newrestype(i,1)=6;
        atomcount=atomcount+1;
        serinfo(sercount,atomcount,1)=x(i,1);
        serinfo(sercount,atomcount,2)=y(i,1);
        serinfo(sercount,atomcount,3)=z(i,1);
        serinfo(sercount,atomcount,4)=105.1;
        if (i==ds || resnumber(i,1)~=resnumber(i+1,1))==1
            sercount=sercount+1;
            rescount=rescount+1;
            res_num(i)=resnumber(i,1);
            resatomcount(rescount,1)=atomcount;
            atomcount=0;
        end
    end
    if strcmp(thr,restype{i})==1
        newrestype(i,1)=7;
        atomcount=atomcount+1;
        thrinfo(thrcount,atomcount,1)=x(i,1);
        thrinfo(thrcount,atomcount,2)=y(i,1);
        thrinfo(thrcount,atomcount,3)=z(i,1);
        thrinfo(thrcount,atomcount,4)=119.1;
        if (i==ds || resnumber(i,1)~=resnumber(i+1,1))==1
            thrcount=thrcount+1;
            rescount=rescount+1;
            res_num(i)=resnumber(i,1);
            resatomcount(rescount,1)=atomcount;
            atomcount=0;
        end
    end
    if strcmp(asp,restype{i})==1
        newrestype(i,1)=8;
        atomcount=atomcount+1;
        aspinfo(aspcount,atomcount,1)=x(i,1);
        aspinfo(aspcount,atomcount,2)=y(i,1);
        aspinfo(aspcount,atomcount,3)=z(i,1);
        aspinfo(aspcount,atomcount,4)=133.1;
        if (i==ds || resnumber(i,1)~=resnumber(i+1,1))==1
            aspcount=aspcount+1;
            rescount=rescount+1;
            res_num(i)=resnumber(i,1);
            resatomcount(rescount,1)=atomcount;
            atomcount=0;
        end
    end
    if strcmp(asn,restype{i})==1 || strcmp('ASH',restype{i})==1
        newrestype(i,1)=9;
        atomcount=atomcount+1;
        asninfo(asncount,atomcount,1)=x(i,1);
        asninfo(asncount,atomcount,2)=y(i,1);
        asninfo(asncount,atomcount,3)=z(i,1);
        asninfo(asncount,atomcount,4)=132.1;
        if (i==ds || resnumber(i,1)~=resnumber(i+1,1))==1
            asncount=asncount+1;
            rescount=rescount+1;
            res_num(i)=resnumber(i,1);
            resatomcount(rescount,1)=atomcount;
            atomcount=0;
        end
    end
    if strcmp(glu,restype{i})==1 || strcmp('GLH',restype{i})==1
        newrestype(i,1)=10;
        atomcount=atomcount+1;
        gluinfo(glucount,atomcount,1)=x(i,1);
        gluinfo(glucount,atomcount,2)=y(i,1);
        gluinfo(glucount,atomcount,3)=z(i,1);
        gluinfo(glucount,atomcount,4)=147.1;
        if (i==ds || resnumber(i,1)~=resnumber(i+1,1))==1
            glucount=glucount+1;
            rescount=rescount+1;
            res_num(i)=resnumber(i,1);
            resatomcount(rescount,1)=atomcount;
            atomcount=0;
        end
    end
    if strcmp(gln,restype{i})==1
        newrestype(i,1)=11;
        atomcount=atomcount+1;
        glninfo(glncount,atomcount,1)=x(i,1);
        glninfo(glncount,atomcount,2)=y(i,1);
        glninfo(glncount,atomcount,3)=z(i,1);
        glninfo(glncount,atomcount,4)=146;
        if (i==ds || resnumber(i,1)~=resnumber(i+1,1))==1
            glncount=glncount+1;
            rescount=rescount+1;
            res_num(i)=resnumber(i,1);
            resatomcount(rescount,1)=atomcount;
            atomcount=0;
        end
    end
    if strcmp(lys,restype{i})==1 || strcmp('LYN',restype{i})==1
        newrestype(i,1)=12;
        atomcount=atomcount+1;
        lysinfo(lyscount,atomcount,1)=x(i,1);
        lysinfo(lyscount,atomcount,2)=y(i,1);
        lysinfo(lyscount,atomcount,3)=z(i,1);
        lysinfo(lyscount,atomcount,4)=146.2;
        if (i==ds || resnumber(i,1)~=resnumber(i+1,1))==1
            lyscount=lyscount+1;
            rescount=rescount+1;
            res_num(i)=resnumber(i,1);
            resatomcount(rescount,1)=atomcount;
            atomcount=0;
        end
    end
    if strcmp(arg,restype{i})==1
        newrestype(i,1)=13;
        atomcount=atomcount+1;
        arginfo(argcount,atomcount,1)=x(i,1);
        arginfo(argcount,atomcount,2)=y(i,1);
        arginfo(argcount,atomcount,3)=z(i,1);
        arginfo(argcount,atomcount,4)=174.2;
        if (i==ds || resnumber(i,1)~=resnumber(i+1,1))==1
            argcount=argcount+1;
            rescount=rescount+1;
            res_num(i)=resnumber(i,1);
            resatomcount(rescount,1)=atomcount;
            atomcount=0;
        end
    end
    if strcmp(cys,restype{i})==1 || strcmp('CYM',restype{i})==1 || strcmp('CYX',restype{i})==1
        newrestype(i,1)=14;
        atomcount=atomcount+1;
        cysinfo(cyscount,atomcount,1)=x(i,1);
        cysinfo(cyscount,atomcount,2)=y(i,1);
        cysinfo(cyscount,atomcount,3)=z(i,1);
        cysinfo(cyscount,atomcount,4)=121.2;
        if (i==ds || resnumber(i,1)~=resnumber(i+1,1))==1
            cyscount=cyscount+1;
            rescount=rescount+1;
            res_num(i)=resnumber(i,1);
            resatomcount(rescount,1)=atomcount;
            atomcount=0;
        end
    end
    if strcmp(met,restype{i})==1
        newrestype(i,1)=15;
        atomcount=atomcount+1;
        metinfo(metcount,atomcount,1)=x(i,1);
        metinfo(metcount,atomcount,2)=y(i,1);
        metinfo(metcount,atomcount,3)=z(i,1);
        metinfo(metcount,atomcount,4)=149.2;
        if (i==ds || resnumber(i,1)~=resnumber(i+1,1))==1
            metcount=metcount+1;
            rescount=rescount+1;
            res_num(i)=resnumber(i,1);
            resatomcount(rescount,1)=atomcount;
            atomcount=0;
        end
    end
    if strcmp(phe,restype{i})==1
        newrestype(i,1)=16;
        atomcount=atomcount+1;
        pheinfo(phecount,atomcount,1)=x(i,1);
        pheinfo(phecount,atomcount,2)=y(i,1);
        pheinfo(phecount,atomcount,3)=z(i,1);
        pheinfo(phecount,atomcount,4)=165.2;
        if (i==ds || resnumber(i,1)~=resnumber(i+1,1))==1
            phecount=phecount+1;
            rescount=rescount+1;
            res_num(i)=resnumber(i,1);
            resatomcount(rescount,1)=atomcount;
            atomcount=0;
        end
    end
    if strcmp(tyr,restype{i})==1
        newrestype(i,1)=17;
        atomcount=atomcount+1;
        tyrinfo(tyrcount,atomcount,1)=x(i,1);
        tyrinfo(tyrcount,atomcount,2)=y(i,1);
        tyrinfo(tyrcount,atomcount,3)=z(i,1);
        tyrinfo(tyrcount,atomcount,4)=181.2;
        if (i==ds || resnumber(i,1)~=resnumber(i+1,1))==1
            tyrcount=tyrcount+1;
            rescount=rescount+1;
            res_num(i)=resnumber(i,1);
            resatomcount(rescount,1)=atomcount;
            atomcount=0;
        end
    end
    if strcmp(trp,restype{i})==1
        newrestype(i,1)=18;
        atomcount=atomcount+1;
        trpinfo(trpcount,atomcount,1)=x(i,1);
        trpinfo(trpcount,atomcount,2)=y(i,1);
        trpinfo(trpcount,atomcount,3)=z(i,1);
        trpinfo(trpcount,atomcount,4)=204.2;
        if (i==ds || resnumber(i,1)~=resnumber(i+1,1))==1
            trpcount=trpcount+1;
            rescount=rescount+1;
            res_num(i)=resnumber(i,1);
            resatomcount(rescount,1)=atomcount;
            atomcount=0;
        end
    end
    if strcmp(his,restype{i})==1 || strcmp('HID',restype{i})==1 || strcmp('HIE',restype{i})==1 || strcmp('HIP',restype{i})==1
        newrestype(i,1)=19;
        atomcount=atomcount+1;
        hisinfo(hiscount,atomcount,1)=x(i,1);
        hisinfo(hiscount,atomcount,2)=y(i,1);
        hisinfo(hiscount,atomcount,3)=z(i,1);
        hisinfo(hiscount,atomcount,4)=155.2;
        if i==ds ||  resnumber(i,1)~=resnumber(i+1,1)
            hiscount=hiscount+1;
            rescount=rescount+1;
            res_num(i)=resnumber(i,1);
            resatomcount(rescount,1)=atomcount;
            atomcount=0;
        end
    end
    if strcmp(pro,restype{i})==1
        newrestype(i,1)=20;
        atomcount=atomcount+1;
        proinfo(procount,atomcount,1)=x(i,1);
        proinfo(procount,atomcount,2)=y(i,1);
        proinfo(procount,atomcount,3)=z(i,1);
        proinfo(procount,atomcount,4)=115.1;
        if i==ds ||  resnumber(i,1)~=resnumber(i+1,1)
            procount=procount+1;
            rescount=rescount+1;
            res_num(i)=resnumber(i,1);
            resatomcount(rescount,1)=atomcount;
            atomcount=0;
        end
    end
    % NEW: only write Cα rows for standard amino acids (newrestype ~= 0)
     if strcmpi(atomtype{i}, 'CA') && newrestype(i) ~= 0
            output(n,1) = newrestype(i);      % residue code 1..20
            output(n,2) = resnumber(i,1);
            output(n,3) = newchainid(i,1);
            output(n,4) = 5;
            output(n,5) = x(i,1);
            output(n,6) = y(i,1);
            output(n,7) = z(i,1);
            n = n + 1;
     end

end
%%
glycount=0;
alacount=0;
valcount=0;
ilecount=0;
leucount=0;
sercount=0;
thrcount=0;
aspcount=0;
asncount=0;
glucount=0;
glncount=0;
lyscount=0;
argcount=0;
cyscount=0;
metcount=0;
phecount=0;
tyrcount=0;
trpcount=0;
hiscount=0;
procount=0;
m=size(output);
totresno=m(1,1);

k=1;
l=1;
err_index=[];
erro=[];
for i=1:totresno-1
    if output(i,1)==1 
		if resatomcount(i,1)>1
			nonerr_index(k)=i;
			k=k+1;
		else
			err_index(l)=i;
            erro(l,:)=[output(i,2) output(i,3)];
			l=l+1;
		end
	end
    if output(i,1)==2 
        if resatomcount(i,1)>4
			nonerr_index(k)=i;
			k=k+1;
		else
			err_index(l)=i;
			erro(l,:)=[output(i,2) output(i,3)];
			l=l+1;
		end
	end
    if output(i,1)==3
        if resatomcount(i,1)>6
			nonerr_index(k)=i;
			k=k+1;
		else
			err_index(l)=i;
erro(l,:)=[output(i,2) output(i,3)];
			l=l+1;
		end
	end
    if output(i,1)==4
        if resatomcount(i,1)>7
			nonerr_index(k)=i;
			k=k+1;
		else
			err_index(l)=i;
erro(l,:)=[output(i,2) output(i,3)];
			l=l+1;
		end
	end
    if output(i,1)==5
        if resatomcount(i,1)>7
			nonerr_index(k)=i;
			k=k+1;
		else
			err_index(l)=i;
erro(l,:)=[output(i,2) output(i,3)];
			l=l+1;
		end
	end
    if output(i,1)==6
        if resatomcount(i,1)>5
			nonerr_index(k)=i;
			k=k+1;
		else
			err_index(l)=i;
erro(l,:)=[output(i,2) output(i,3)];
			l=l+1;
		end
	end
    if output(i,1)==7
        if resatomcount(i,1)>5
			nonerr_index(k)=i;
			k=k+1;
		else
			err_index(l)=i;
erro(l,:)=[output(i,2) output(i,3)];
			l=l+1;
		end
	end
    if output(i,1)==8
        if resatomcount(i,1)>7
			nonerr_index(k)=i;
			k=k+1;
		else
			err_index(l)=i;
erro(l,:)=[output(i,2) output(i,3)];
			l=l+1;
		end
	end
    if output(i,1)==9
        if resatomcount(i,1)>7
			nonerr_index(k)=i;
			k=k+1;
		else
			err_index(l)=i;
erro(l,:)=[output(i,2) output(i,3)];
			l=l+1;
		end
	end
    if output(i,1)==10
        if resatomcount(i,1)>8
			nonerr_index(k)=i;
			k=k+1;
		else
			err_index(l)=i;
erro(l,:)=[output(i,2) output(i,3)];
			l=l+1;
		end
	end
    if output(i,1)==11
        if resatomcount(i,1)>8
			nonerr_index(k)=i;
			k=k+1;
		else
			err_index(l)=i;
erro(l,:)=[output(i,2) output(i,3)];
			l=l+1;
		end
	end
    if output(i,1)==12
        if resatomcount(i,1)>8
			nonerr_index(k)=i;
			k=k+1;
		else
			err_index(l)=i;
erro(l,:)=[output(i,2) output(i,3)];
			l=l+1;
		end
	end
    if output(i,1)==13
        if resatomcount(i,1)>10
			nonerr_index(k)=i;
			k=k+1;
		else
			err_index(l)=i;
erro(l,:)=[output(i,2) output(i,3)];
			l=l+1;
		end
	end
    if output(i,1)==14
        if resatomcount(i,1)>5
			nonerr_index(k)=i;
			k=k+1;
		else
			err_index(l)=i;
erro(l,:)=[output(i,2) output(i,3)];
			l=l+1;
		end
	end
    if output(i,1)==15
        if resatomcount(i,1)>6
			nonerr_index(k)=i;
			k=k+1;
		else
			err_index(l)=i;
erro(l,:)=[output(i,2) output(i,3)];
			l=l+1;
		end
	end
    if output(i,1)==16
        if resatomcount(i,1)>10
			nonerr_index(k)=i;
			k=k+1;
		else
			err_index(l)=i;
erro(l,:)=[output(i,2) output(i,3)];
			l=l+1;
		end
	end
    if output(i,1)==17
        if resatomcount(i,1)>11
			nonerr_index(k)=i;
			k=k+1;
		else
			err_index(l)=i;
erro(l,:)=[output(i,2) output(i,3)];
			l=l+1;
		end
	end
    if output(i,1)==18
        if resatomcount(i,1)>13
			nonerr_index(k)=i;
			k=k+1;
		else
			err_index(l)=i;
erro(l,:)=[output(i,2) output(i,3)];
			l=l+1;
		end
	end
    if output(i,1)==19
        if resatomcount(i,1)>9
			nonerr_index(k)=i;
			k=k+1;
		else
			err_index(l)=i;
erro(l,:)=[output(i,2) output(i,3)];
			l=l+1;
		end
	end
    if output(i,1)==20
        if resatomcount(i,1)>6
			nonerr_index(k)=i;
			k=k+1;
		else
			err_index(l)=i;
erro(l,:)=[output(i,2) output(i,3)];
			l=l+1;
		end
	end
end
%%
for i=1:totresno
	if output(i,1)==1
        if resatomcount(i,1)>1
        glycount=glycount+1;
        output(i,8)=glyinfo(glycount,2,1);
        output(i,9)=glyinfo(glycount,2,2);
        output(i,10)=glyinfo(glycount,2,3);
        output(i,11)=glyinfo(glycount,1,4);
        end
    end
    if output(i,1)==2
        if resatomcount(i,1)>4
        alacount=alacount+1;
        output(i,8)=alainfo(alacount,5,1);
        output(i,9)=alainfo(alacount,5,2);
        output(i,10)=alainfo(alacount,5,3);
        output(i,11)=alainfo(alacount,1,4);
        else
            alacount=alacount+1;
        output(i,8)=output(i,5);
        output(i,9)=output(i,6);
        output(i,10)=output(i,7);
        output(i,11)=alainfo(alacount,1,4);    
        end
    end
    if output(i,1)==3
        if resatomcount(i,1)>6
        valcount=valcount+1;
        output(i,8)=(valinfo(valcount,6,1)+valinfo(valcount,7,1))/2;
        output(i,9)=(valinfo(valcount,6,2)+valinfo(valcount,7,2))/2;
        output(i,10)=(valinfo(valcount,6,3)+valinfo(valcount,7,3))/2;
        output(i,11)=valinfo(valcount,1,4);
        else
            valcount=valcount+1;
        output(i,8)=output(i,5);
        output(i,9)=output(i,6);
        output(i,10)=output(i,7);
        output(i,11)=valinfo(valcount,1,4); 
        end
    end
    if output(i,1)==4
        if resatomcount(i,1)>7
        ilecount=ilecount+1;
        output(i,8)=ileinfo(ilecount,8,1);
        output(i,9)=ileinfo(ilecount,8,2);
        output(i,10)=ileinfo(ilecount,8,3);
        output(i,11)=ileinfo(ilecount,1,4);
        else
            ilecount=ilecount+1;
        output(i,8)=output(i,5);
        output(i,9)=output(i,6);
        output(i,10)=output(i,7);
        output(i,11)=ileinfo(ilecount,1,4); 
        end
    end
    if output(i,1)==5
        if resatomcount(i,1)>7
        leucount=leucount+1;
        output(i,8)=(leuinfo(leucount,7,1)+leuinfo(leucount,8,1))/2;
        output(i,9)=(leuinfo(leucount,7,2)+leuinfo(leucount,8,2))/2;
        output(i,10)=(leuinfo(leucount,7,3)+leuinfo(leucount,8,3))/2;
        output(i,11)=leuinfo(leucount,1,4);
        else
            leucount=leucount+1;
        output(i,8)=output(i,5);
        output(i,9)=output(i,6);
        output(i,10)=output(i,7);
        output(i,11)=leuinfo(leucount,1,4); 
        end
    end
    if output(i,1)==6
        if resatomcount(i,1)>5
        sercount=sercount+1;
        output(i,8)=serinfo(sercount,6,1);
        output(i,9)=serinfo(sercount,6,2);
        output(i,10)=serinfo(sercount,6,3);
        output(i,11)=serinfo(sercount,1,4);
        else
            sercount=sercount+1;
        output(i,8)=output(i,5);
        output(i,9)=output(i,6);
        output(i,10)=output(i,7);
        output(i,11)=serinfo(sercount,1,4); 
        end
    end
    if output(i,1)==7
        if resatomcount(i,1)>5
        thrcount=thrcount+1;
        output(i,8)=thrinfo(thrcount,6,1);
        output(i,9)=thrinfo(thrcount,6,2);
        output(i,10)=thrinfo(thrcount,6,3);
        output(i,11)=thrinfo(thrcount,1,4);
        else
            thrcount=thrcount+1;
        output(i,8)=output(i,5);
        output(i,9)=output(i,6);
        output(i,10)=output(i,7);
        output(i,11)=thrinfo(thrcount,1,4); 
        end
    end
    if output(i,1)==8
        if resatomcount(i,1)>7
        aspcount=aspcount+1;
        output(i,8)=(aspinfo(aspcount,7,1)+aspinfo(aspcount,8,1))/2;
        output(i,9)=(aspinfo(aspcount,7,2)+aspinfo(aspcount,8,2))/2;
        output(i,10)=(aspinfo(aspcount,7,3)+aspinfo(aspcount,8,3))/2;
        output(i,11)=aspinfo(aspcount,1,4);
        else
            aspcount=aspcount+1;
        output(i,8)=output(i,5);
        output(i,9)=output(i,6);
        output(i,10)=output(i,7);
        output(i,11)=aspinfo(aspcount,1,4); 
        end
    end
    if output(i,1)==9
        if resatomcount(i,1)>7
        asncount=asncount+1;
        output(i,8)=(asninfo(asncount,7,1)+asninfo(asncount,8,1))/2;
        output(i,9)=(asninfo(asncount,7,2)+asninfo(asncount,8,2))/2;
        output(i,10)=(asninfo(asncount,7,3)+asninfo(asncount,8,3))/2;
        output(i,11)=asninfo(asncount,1,4);
        else
            asncount=asncount+1;
        output(i,8)=output(i,5);
        output(i,9)=output(i,6);
        output(i,10)=output(i,7);
        output(i,11)=asninfo(asncount,1,4); 
        end
    end
    if output(i,1)==10
        if resatomcount(i,1)>8
        glucount=glucount+1;
        output(i,8)=(gluinfo(glucount,8,1)+gluinfo(glucount,9,1))/2;
        output(i,9)=(gluinfo(glucount,8,2)+gluinfo(glucount,9,2))/2;
        output(i,10)=(gluinfo(glucount,8,3)+gluinfo(glucount,9,3))/2;
        output(i,11)=gluinfo(glucount,1,4);
        else
            glucount=glucount+1;
        output(i,8)=output(i,5);
        output(i,9)=output(i,6);
        output(i,10)=output(i,7);
        output(i,11)=gluinfo(glucount,1,4); 
        end
    end
    if output(i,1)==11
        if resatomcount(i,1)>8
        glncount=glncount+1;
        output(i,8)=(glninfo(glncount,8,1)+glninfo(glncount,9,1))/2;
        output(i,9)=(glninfo(glncount,8,2)+glninfo(glncount,9,2))/2;
        output(i,10)=(glninfo(glncount,8,3)+glninfo(glncount,9,3))/2;
        output(i,11)=glninfo(glncount,1,4);
        else
            glncount=glncount+1;
        output(i,8)=output(i,5);
        output(i,9)=output(i,6);
        output(i,10)=output(i,7);
        output(i,11)=glninfo(glncount,1,4); 
        end
    end
    if output(i,1)==12
        if resatomcount(i,1)>8
        lyscount=lyscount+1;
        output(i,8)=lysinfo(lyscount,9,1);
        output(i,9)=lysinfo(lyscount,9,2);
        output(i,10)=lysinfo(lyscount,9,3);
        output(i,11)=lysinfo(lyscount,1,4);
        else
            lyscount=lyscount+1;
        output(i,8)=output(i,5);
        output(i,9)=output(i,6);
        output(i,10)=output(i,7);
        output(i,11)=lysinfo(lyscount,1,4); 
        end
    end
    if output(i,1)==13
        if resatomcount(i,1)>10
        argcount=argcount+1;
        output(i,8)=(arginfo(argcount,8,1)+arginfo(argcount,10,1)+arginfo(argcount,11,1))/3;
        output(i,9)=(arginfo(argcount,8,2)+arginfo(argcount,10,2)+arginfo(argcount,11,2))/3;
        output(i,10)=(arginfo(argcount,8,3)+arginfo(argcount,10,3)+arginfo(argcount,11,3))/3;
        output(i,11)=arginfo(argcount,1,4);
        else
            argcount=argcount+1;
        output(i,8)=output(i,5);
        output(i,9)=output(i,6);
        output(i,10)=output(i,7);
        output(i,11)=arginfo(argcount,1,4); 
        end
    end
    if output(i,1)==14
        if resatomcount(i,1)>5
        cyscount=cyscount+1;
        output(i,8)=cysinfo(cyscount,6,1);
        output(i,9)=cysinfo(cyscount,6,2);
        output(i,10)=cysinfo(cyscount,6,3);
        output(i,11)=cysinfo(cyscount,1,4);
        else
            cyscount=cyscount+1;
        output(i,8)=output(i,5);
        output(i,9)=output(i,6);
        output(i,10)=output(i,7);
        output(i,11)=cysinfo(cyscount,1,4); 
        end
    end
    if output(i,1)==15
        if resatomcount(i,1)>6
        metcount=metcount+1;
        output(i,8)=metinfo(metcount,7,1);
        output(i,9)=metinfo(metcount,7,2);
        output(i,10)=metinfo(metcount,7,3);
        output(i,11)=metinfo(metcount,1,4);
        else
            metcount=metcount+1;
        output(i,8)=output(i,5);
        output(i,9)=output(i,6);
        output(i,10)=output(i,7);
        output(i,11)=metinfo(metcount,1,4); 
        end
    end
    if output(i,1)==16
        if resatomcount(i,1)>10
        phecount=phecount+1;
        output(i,8)=(pheinfo(phecount,6,1)+pheinfo(phecount,7,1)+pheinfo(phecount,8,1)+pheinfo(phecount,9,1)+pheinfo(phecount,10,1)+pheinfo(phecount,11,1))/6;
        output(i,9)=(pheinfo(phecount,6,2)+pheinfo(phecount,7,2)+pheinfo(phecount,8,2)+pheinfo(phecount,9,2)+pheinfo(phecount,10,2)+pheinfo(phecount,11,2))/6;
        output(i,10)=(pheinfo(phecount,6,3)+pheinfo(phecount,7,3)+pheinfo(phecount,8,3)+pheinfo(phecount,9,3)+pheinfo(phecount,10,3)+pheinfo(phecount,11,3))/6;
        output(i,11)=pheinfo(phecount,1,4);
        else
            phecount=phecount+1;
        output(i,8)=output(i,5);
        output(i,9)=output(i,6);
        output(i,10)=output(i,7);
        output(i,11)=pheinfo(phecount,1,4); 
        end
    end
    if output(i,1)==17
        if resatomcount(i,1)>11
        tyrcount=tyrcount+1;
        output(i,8)=(tyrinfo(tyrcount,6,1)+tyrinfo(tyrcount,7,1)+tyrinfo(tyrcount,8,1)+tyrinfo(tyrcount,9,1)+tyrinfo(tyrcount,10,1)+tyrinfo(tyrcount,11,1)+tyrinfo(tyrcount,12,1))/7;
        output(i,9)=(tyrinfo(tyrcount,6,2)+tyrinfo(tyrcount,7,2)+tyrinfo(tyrcount,8,2)+tyrinfo(tyrcount,9,2)+tyrinfo(tyrcount,10,2)+tyrinfo(tyrcount,11,2)+tyrinfo(tyrcount,12,2))/7;
        output(i,10)=(tyrinfo(tyrcount,6,3)+tyrinfo(tyrcount,7,3)+tyrinfo(tyrcount,8,3)+tyrinfo(tyrcount,9,3)+tyrinfo(tyrcount,10,3)+tyrinfo(tyrcount,11,3)+tyrinfo(tyrcount,12,3))/7;
        output(i,11)=tyrinfo(tyrcount,1,4);
        else
            tyrcount=tyrcount+1;
        output(i,8)=output(i,5);
        output(i,9)=output(i,6);
        output(i,10)=output(i,7);
        output(i,11)=tyrinfo(tyrcount,1,4); 
        end
    end
    if output(i,1)==18
        if resatomcount(i,1)>13
        trpcount=trpcount+1;
        output(i,8)=(trpinfo(trpcount,7,1)+trpinfo(trpcount,8,1)+trpinfo(trpcount,9,1)+trpinfo(trpcount,10,1)+trpinfo(trpcount,11,1)+trpinfo(trpcount,12,1)+trpinfo(trpcount,13,1)+trpinfo(trpcount,14,1))/8;
        output(i,9)=(trpinfo(trpcount,7,2)+trpinfo(trpcount,8,2)+trpinfo(trpcount,9,2)+trpinfo(trpcount,10,2)+trpinfo(trpcount,11,2)+trpinfo(trpcount,12,2)+trpinfo(trpcount,13,2)+trpinfo(trpcount,14,2))/8;
        output(i,10)=(trpinfo(trpcount,7,3)+trpinfo(trpcount,8,3)+trpinfo(trpcount,9,3)+trpinfo(trpcount,10,3)+trpinfo(trpcount,11,3)+trpinfo(trpcount,12,3)+trpinfo(trpcount,13,3)+trpinfo(trpcount,14,3))/8;
        output(i,11)=trpinfo(trpcount,1,4);
        else
            trpcount=trpcount+1;
        output(i,8)=output(i,5);
        output(i,9)=output(i,6);
        output(i,10)=output(i,7);
        output(i,11)=trpinfo(trpcount,1,4); 
        end
    end
    if output(i,1)==19
        if resatomcount(i,1)>9
        hiscount=hiscount+1;
        output(i,8)=(hisinfo(hiscount,6,1)+hisinfo(hiscount,7,1)+hisinfo(hiscount,8,1)+hisinfo(hiscount,9,1)+hisinfo(hiscount,10,1))/5;
        output(i,9)=(hisinfo(hiscount,6,2)+hisinfo(hiscount,7,2)+hisinfo(hiscount,8,2)+hisinfo(hiscount,9,2)+hisinfo(hiscount,10,2))/5;
        output(i,10)=(hisinfo(hiscount,6,3)+hisinfo(hiscount,7,3)+hisinfo(hiscount,8,3)+hisinfo(hiscount,9,3)+hisinfo(hiscount,10,3))/5;
        output(i,11)=hisinfo(hiscount,1,4);
        else
            hiscount=hiscount+1;
        output(i,11)=hisinfo(hiscount,1,4);
        output(i,8)=output(i,5);
        output(i,9)=output(i,6);
        output(i,10)=output(i,7);
        output(i,11)=hisinfo(hiscount,1,4); 
        end
    end
    if output(i,1)==20
        if resatomcount(i,1)>6
        procount=procount+1;
        output(i,8)=(proinfo(procount,5,1)+proinfo(procount,6,1)+proinfo(procount,7,1))/3;
        output(i,9)=(proinfo(procount,5,2)+proinfo(procount,6,2)+proinfo(procount,7,2))/3;
        output(i,10)=(proinfo(procount,5,3)+proinfo(procount,6,3)+proinfo(procount,7,3))/3;
        output(i,11)=proinfo(procount,1,4);
        else
            procount=procount+1;
        output(i,8)=output(i,5);
        output(i,9)=output(i,6);
        output(i,10)=output(i,7);
        output(i,11)=proinfo(procount,1,4); 
        end
    end
end
%%
h=size(output);
resno=h(1);
k=0;
b=0;
for i=1:ds-1
    b=b+1;
    if strcmp('N',atomtype{i+1})~=0 || strcmp(restype{i},restype{i+1})==0 %|| strcmp('OXT',atomtype(i+1,:))==1
        k=k+1;
        output(k,4)=b;
        b=0;
    end
end
output(k+1,4)=b+1;


for i=1:resno
    result(i,1)=output(i,2);
    result(i,2)=output(i,1);
    result(i,3)=output(i,5);
    result(i,4)=output(i,6);
    result(i,5)=output(i,7);
    result(i,6)=output(i,8);
    result(i,7)=output(i,9);
    result(i,8)=output(i,10);
    result(i,9)=output(i,4);
    result(i,10)=output(i,3);
    result(i,11)=output(i,11);
end




%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lacking chain cooridnates are replaced by backbone coordinates
for i=1:resno
    if result(i,6)==0 & result(i,7)==0 & result(i,8)==0
        result(i,6)=result(i,3);
        result(i,7)=result(i,4);
        result(i,8)=result(i,5);
    end
end

    % === Keep only the selected chain (by rank) ===
if ~isempty(chainID)
    if isempty(rank)
        warning('Requested chainID "%s" not found. No filtering applied.', chainID);
    else
        keep = (result(:,10) == rank);
        result = result(keep, :);
        % (Optional) normalize chain column to 1 after filtering:
        % result(:,10) = 1;
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BETA-CARBON COORDINATES (for use in Thomas-Dill)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ATOMISTIC COORDINATES


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ANGLE CALCULATION



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WRITING OUTPUT FILES
File=strcat( ProteinName , '.cor' );
%File = ['../runs/' ProteinName '.cor'];
fid=fopen(File,'w');
fprintf(fid, '%5.0f %4.0f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %4.0f %2.0f %4.1f\n', result');
fclose(fid);
if length(erro)>0
File = 'problem';
fid = fopen(File,'w');
fprintf(fid,'There are problems about number of atoms on following residues\n');
for i=1:size(erro,1)
    fprintf(fid,'%d of %d chain\n',erro(i,1),erro(i,2));
end
fclose(fid);
end
end