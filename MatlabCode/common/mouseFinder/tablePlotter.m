function tablePlotter

num = evalin('base','num');raw = evalin('base','raw');txt = evalin('base','txt');

% columns in excel file
ID = 'A';
Genotype = 'H';
User = 'I';
ProjectI_V = 'I';
DOB = 'C';
Sex = 'B';


% specify min and max age
ages = getAges(raw,DOB);
sexes = getSelectedCol(raw,Sex);
geno = getSelectedCol(raw,Genotype);
users = getSelectedCol(raw,User);

selections = {'user','Sam'};
selections = {'age',[0 4]};

% selections = {'age',[0 6]};
% selections = {'user','Aras','age',[3 6]};

ind = 1;
for ii = 1:2:length(selections)
    selType = selections{ii};
    selArg = selections{ii+1};
    if strcmp(selType,'user')
        IndexC = strfind(users,selArg);
        selection{ind} = find(not(cellfun('isempty', IndexC)));
    end
    if strcmp(selType,'age')
        selection{ind} = find(ages < selArg(2) & ages > selArg(1));
    end
    if strcmp(selType,'geno')
        IndexC = strfind(geno,selArg);
        selection{ind} = find(not(cellfun('isempty', IndexC)));
    end
    if strcmp(selType,'sex')
        selection{ind} = strfind(sexes',selArg);
    end
    ind = ind + 1;
end

if exist('selection','var')
    common = selection{1};
    if length(selection) > 1
        for ii = 2:length(selection)
            common = intersect(common,selection{ii});
        end
    end
else
    common = 1:size(users,1);
end

ids = raw(:,1);
varNames = {'No','ID','Gender','Age','User','GenoType'};
No = (1:length(common))';
T = table(No,ids(common),sexes(common),ages(common)',users(common),geno(common),'VariableNames',varNames)
temp = sexes(common);
totalAnimals = length(common);
males = length(strfind(temp','M'));
females = totalAnimals - males;
display(sprintf('Total number of Animals = %d, (%dM, %dF)',totalAnimals,males,females));


function ages = getAges(raw,colA,selRows)
Alphabets = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
col = strfind(Alphabets,colA);
if ~exist('selRows','var')
    selRows = 1:size(raw,1);
end
for ii = 1:length(selRows)
    if ii == 103
        n = 0;
    end
    thisAge = raw{selRows(ii),col};
    numdays = today - datenum(thisAge);
    ages(ii) = numdays/30;
%         [Y,M,D,H,MN,S] = datevec(numdays);
%         ages(ii) = (Y * 12) + (M-1) + (D/30);
end


function rows = getSelectedCol(raw,colA,selRows)
Alphabets = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
col = strfind(Alphabets,colA);
if ~exist('selRows','var')
    selRows = 1:size(raw,1);
end
for ii = 1:length(selRows)
    if isnan(raw{selRows(ii),col})
        if strcmp(colA,'B')
            rows(ii,1) = ' ';
        else
            rows{ii,1} = ' ';
        end
    else
        if strcmp(colA,'B')
            rows(ii,1) = raw{selRows(ii),col};
        else
            rows{ii,1} = raw{selRows(ii),col};
        end
    end
end
