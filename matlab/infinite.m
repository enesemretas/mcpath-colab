%%%%%%%%%%%%%%%%INFINITE PATH GENERATOR%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generation of a long path with a huge step number entered by the user

function y = infinite(ProteinName, path_length)
y = 0;

% --- select potential type (kept as in original) ---
pottype = '1';
if pottype=='1'
    Type = '_atomistic'; tip = 'atomistic';
elseif pottype=='2'
    Type = '_BJ';        tip = 'BJ';
elseif pottype=='3'
    Type = '_TD';        tip = 'TD';
elseif pottype=='4'
    Type = '_dist';      tip = 'dist';
elseif pottype=='5'
    Type = '_markov';    tip = 'markov';
else
    Type = '_atomistic'; tip = 'atomistic';
end

% --- read inputs (.cor and probability .out) ---
File = strcat(ProteinName, '.cor');
if ~exist(File, 'file')
    error('Missing file: %s', File);
end
coor = dlmread(File);  % input, output of the readpdb.f

File = strcat(ProteinName, Type, '.out');
if ~exist(File, 'file')
    error('Missing file: %s', File);
end
normalized = dlmread(File); % probability matrix

% --- starting indices and path length (kept default behavior) ---
baslangic  = 1;  % starting residue number
baslangic1 = 1;  % chain number

% robust parse for Octave and MATLAB (keeps your original intent)
if isnumeric(path_length)
    pathlength = path_length;
else
    pathlength = str2num(path_length); %#ok<ST2NM>
end
if isempty(pathlength) || ~isscalar(pathlength) || pathlength<1
    error('Invalid path_length. Provide a positive scalar or numeric string.');
end
pathlength = floor(pathlength);

fromres = [num2str(baslangic) ',' num2str(baslangic1)]; %#ok<NASGU>

m = size(normalized);
resno = m(1);
upper = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find initial row index from coor (as in original)
initial = 1;
for i = 1:resno   % correcting the input residue numbers
    if coor(i,10)==baslangic1 && coor(i,1)==baslangic  % to obtain row numbers
        initial = i;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build newmtrx with header row/col (as in original)
newmtrx = zeros(resno+1, resno+1);
for i = 1:resno+1
    newmtrx(1,i) = i-1;
    newmtrx(i,1) = i-1;
end
for i = 1:resno
    for j = 1:resno
        newmtrx(i+1,j+1) = normalized(i,j);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NEW MATRIX GENERATION (kept identical; minor safety for sums)
if str2num(pottype)==1 %#ok<ST2NM>
    cap = ones(resno,resno); %#ok<NASGU>
    sumv = 0;
    for i = 2:resno+1
        for j = 2:resno+1
            sumv = newmtrx(i,j) + sumv;
        end
        if sumv==0
            % avoid division by zero; keep row as zeros (no outgoing probs)
            % preserves original behavior as closely as possible
        else
            for j = 2:resno+1
                newmtrx(i,j) = newmtrx(i,j)/sumv;
                if i==j, newmtrx(i,j) = 1; end
            end
        end
        minimum = min(newmtrx(i,2:end));
        for j = 2:resno+1
            if minimum==newmtrx(i,j)
                newmtrx(i,j) = 0;
            end
            if i==j
                newmtrx(i,j) = 0;
            end
        end
        sumv = 0;
        for j = 2:resno+1
            sumv = newmtrx(i,j) + sumv;
        end
        if sumv==0
            % keep row zeros if no transitions remain
        else
            for j = 2:resno+1
                newmtrx(i,j) = newmtrx(i,j)/sumv;
                if i==j, newmtrx(i,j) = 0; end
            end
        end
        sumv = 0;
    end
    lower = 0;
else
    cap = ones(resno,resno); %#ok<NASGU>
    sumv = 0;
    for i = 2:resno+1
        for j = 2:resno+1
            sumv = newmtrx(i,j) + sumv;
        end
        average = sumv / resno;
        if sumv==0
            % leave row zeros
        else
            for j = 2:resno+1
                newmtrx(i,j) = newmtrx(i,j)/sumv;
                if i==j, newmtrx(i,j) = 1; end
            end
        end
        for j = 2:resno+1
            if newmtrx(i,j) < average
                newmtrx(i,j) = 0;
            end
            if i==j
                newmtrx(i,j) = 0;
            end
        end
        sumv = 0;
        for j = 2:resno+1
            sumv = newmtrx(i,j) + sumv;
        end
        counter = 0; %#ok<NASGU>
        if sumv==0
            % keep zeros
        else
            for j = 2:resno+1
                newmtrx(i,j) = newmtrx(i,j)/sumv;
                if i==j
                    newmtrx(i,j) = 0;
                end
                if newmtrx(i,j) ~= 0
                    counter = counter+1;
                    cap(i-1,counter) = newmtrx(i,j); %#ok<AGROW>
                end
            end
        end
        sumv = 0;
    end
    lower = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PATH GENERATION (kept same; just minor safety on empty transitions)
runs = 0;          % counts the number of paths generated
prob = 0;          % probability cost accumulator (unused but kept)
taban = max(1, floor(pathlength/10));
permtrx = newmtrx;
runsonscreen = 0;  %#ok<NASGU>

startres = 1;         % original sets 1
runs = runs+1;
newmtrx = permtrx;
counter = 1;
sequence(counter, runs) = startres; %#ok<AGROW>
path(runs, counter) = coor(startres,1) + coor(startres,10)/10; %#ok<AGROW>

while counter ~= pathlength
    counter = counter + 1;
    j = 0;    % number of nonzero energy residues
    clear gen
    clear P
    clear Q
    P = 0;
    for i = 2:resno+1
        if newmtrx(startres+1, i) > lower && newmtrx(startres+1, i) < upper
            j = j + 1;
            P(1,j) = newmtrx(1,i);
            P(2,j) = newmtrx(startres+1,i);
        end
    end
    P = abs(P);

    sumv = 0;
    for i = 1:j
        sumv = sumv + P(2,i);
    end

    if j==0 || sumv==0
        % No outgoing transitions: stay in place (keeps original spirit without crashing)
        startres = startres;
    else
        for i = 1:j
            Q(1,i) = P(1,i);
            Q(2,i) = P(2,i)/sumv;
        end
        number = rand;
        sumv = 0;
        for i = 1:j
            gen(1,i) = Q(1,i);
            sumv = sumv + Q(2,i);
            gen(2,i) = sumv;
        end
        for i = 1:j
            if number < gen(2,i)
                startres = gen(1,i);
                break;
            end
        end
    end

    sequence(counter, runs) = startres; %#ok<AGROW>
    path(runs, counter) = coor(sequence(counter,runs),1) + coor(sequence(counter,runs),10)/10; %#ok<AGROW>

    if mod(counter, taban) == 0
        fprintf('Number of steps generated so far: %d\n', counter);
    end
end

s = size(path);
pl = s(2); % total path length
pn = s(1); %#ok<NASGU> % total number of paths

for i = 1:pl
    if max(abs(path(:, (pl+1-i)))) ~= 0
        j = pl+1-i; %#ok<NASGU>
        break
    end
end
paths = path; % original "trim" logic preserved but effectively already tight

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WRITING PATHS
pl_str = num2str(pl);
what = strcat(ProteinName, Type, '_', pl_str, 'steps_infinite');
File = [what '.path'];
dlmwrite(File, path, '\t');

pl_num = str2num(pl_str); %#ok<ST2NM>
coor(:,1) = coor(:,1) + coor(:,10)/10;
stats = coor(:,1);
stats(:,2) = 0;
for i = 1:pl_num
    x = find(stats(:,1) == path(1,i));
    if ~isempty(x)
        stats(x,2) = stats(x,2) + 1;
    end
end

end
