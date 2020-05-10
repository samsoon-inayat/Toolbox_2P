function list = selectCells16(selAnimals,type,adata)

% adata = evalin('base','data');
mData = evalin('base','mData');
colors = mData.colors;
sigColor = mData.sigColor;
% selAnimals = 1:3;
% selAnimals = 6:11;
mData.belt_length = adata{selAnimals(1)}{1}{1}.belt_length;
n = 0;

%%
for jj = 1:length(selAnimals)
    for ii = 1:7%length(data)
         [tempD cnsjj] = getVariableValues(adata{selAnimals(jj)},'placeCells5',ii);
         distD{jj,ii} = tempD;
         cns{jj,ii} = cnsjj;
         pcs(jj,ii) = sum(tempD);
    end
    numCells(jj,1) = length(tempD);
end
for jj = 1:length(selAnimals)
    allpcsU = distD{jj,1};  allpcsA = distD{jj,1};
    for ii = 2:7%length(data)
         allpcsU = allpcsU | distD{jj,ii}; allpcsA = allpcsA & distD{jj,ii};
         lastPCs = distD{jj,ii-1}; currentPCs = distD{jj,ii};
         remained(jj,ii) = 100*sum(lastPCs & currentPCs)/sum(lastPCs);
         disrupted(jj,ii) = 100*sum(lastPCs & ~currentPCs)/sum(lastPCs);
         newones(jj,ii) = 100*sum(~lastPCs & currentPCs)/sum(currentPCs);
         remainedC{jj,ii} = (lastPCs & currentPCs);
         disruptedC{jj,ii} = (lastPCs & ~currentPCs);
         newonesC{jj,ii} = (~lastPCs & currentPCs);
    end
    allpcsUall{jj} = allpcsU;
    allpcsAall{jj} = allpcsA;
end

for ii = 1:7%length(data)
    distDi = [];
    for jj = 1:length(selAnimals)
        [tempD cnsjj] = getVariableValues(adata{selAnimals(jj)},'placeCells3',ii);
         distDi = [distDi;tempD'];
    end
    allpcscontext{ii} = distDi;
end

if ~isempty(strfind(type,'Common'))
    if strcmp(type,'Common') % to get cells common in all contexts
        commonCells = allpcscontext{1};
        for ii = 2:7
            commonCells = commonCells & allpcscontext{ii};
        end
        list = commonCells;
    else
        numberOfCs = length(type) - length('Common');
        st = length('Common') + 1;
        fc = str2double(type(st)); 
        commonCells = allpcscontext{fc};
        for ii = (st+1):length(type)
            fc = str2double(type(ii)); 
            commonCells = commonCells & allpcscontext{fc};
        end
        list = commonCells;
    end
    return;
end

if ~isempty(strfind(type,'Only')) % to get cells only in a particular context
    sc = str2double(type(end)); 
    osc = setdiff(1:7,sc);
    commonCells = allpcscontext{sc};
    for ii = 1:length(osc)
        commonCells = commonCells & ~allpcscontext{osc(ii)};
    end
    list = commonCells;
    return
end

if ~isempty(strfind(type,'New')) % to get new cells in a particular context
    sc = str2double(type(end)); 
    lastpcs = allpcscontext{sc-1};
    currentpcs = allpcscontext{sc};
    list = currentpcs & ~lastpcs;
    return
end

if ~isempty(strfind(type,'Disrupted')) % to get disruptedcells in a particular context
    sc = str2double(type(end)); 
    nextpcs = allpcscontext{sc+1};
    currentpcs = allpcscontext{sc};
    list = currentpcs & ~nextpcs;
    return
end

if ~isempty(strfind(type,'Remained')) % to get disruptedcells in a particular context
    sc = str2double(type(end)); 
    nextpcs = allpcscontext{sc+1};
    currentpcs = allpcscontext{sc};
    list = currentpcs & nextpcs;
    return
end

