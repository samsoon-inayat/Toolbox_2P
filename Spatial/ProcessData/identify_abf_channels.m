function channels = identify_abf_channels(d,si)



% db.channel{1} = 'frames';
% db.channel{2} = 'ch_a';
% db.channel{3} = 'ch_b';
% db.channel{4} = 'photo_sensor';
% db.channel{5} = 'air_puff';

if size(d,2) == 4
    channels = cell(1,4);
    for ii = 1:size(d,2)
        thisSignal = d(:,ii);
        risingEdges = find_rising_edge(d(:,ii),0.5,2);
        fallingEdges = find_falling_edge(d(:,ii),-0.5,2);
        edges(ii,1) = length(risingEdges);
        edges(ii,2) = length(fallingEdges);
        vals(ii,1) = sum(thisSignal > 2.5);
        vals(ii,2) = sum(thisSignal < 2.5);
        ratioV(ii) = vals(ii,1)/vals(ii,2);
    end
    
    channels{1} = 'frames';

    ind = find(ratioV(2:end) == max(ratioV(2:end)));
    channels{ind+1} = 'photo_sensor';

    inds = setdiff([1:size(d,2)],[1 (ind+1)]);

    ind = min(inds);
    channels{ind} = 'ch_a';
    ind = setdiff(inds,ind);
    channels{ind} = 'ch_b';
end


if size(d,2) == 5
    channels = cell(1,5);
    for ii = 1:size(d,2)
        thisSignal = d(:,ii);
        risingEdges = find_rising_edge(d(:,ii),0.5,2);
        fallingEdges = find_falling_edge(d(:,ii),-0.5,2);
        edges(ii,1) = length(risingEdges);
        edges(ii,2) = length(fallingEdges);
        vals(ii,1) = sum(thisSignal > 2.5);
        vals(ii,2) = sum(thisSignal < 2.5);
        ratioV(ii) = vals(ii,1)/vals(ii,2);
    end
    
    channels{1} = 'frames';

    ind = find(ratioV(2:end) == max(ratioV(2:end)));
    channels{ind+1} = 'photo_sensor';

    inds = setdiff([1:size(d,2)],[1 (ind+1)]);

    ind = find(edges(inds,1) == min(edges(inds,1)));

    channels{inds(ind)} = 'air_puff';

    inds = setdiff(inds,inds(ind));

    ind = min(inds);
    channels{ind} = 'ch_a';
    ind = setdiff(inds,ind);
    channels{ind} = 'ch_b';
end

if size(d,2) == 6
    channels = cell(1,6);
    for ii = 1:size(d,2)
        thisSignal = d(:,ii);
        risingEdges = find_rising_edge(d(:,ii),0.5,2);
        fallingEdges = find_falling_edge(d(:,ii),-0.5,2);
        edges(ii,1) = length(risingEdges);
        edges(ii,2) = length(fallingEdges);
        vals(ii,1) = sum(thisSignal > 2.5);
        vals(ii,2) = sum(thisSignal < 2.5);
        ratioV(ii) = vals(ii,1)/vals(ii,2);
    end
    
    channels{1} = 'frames';

    ind = find(ratioV(2:end) == max(ratioV(2:end)));
    channels{ind+1} = 'photo_sensor';
    inds = setdiff([1:size(d,2)],[1 (ind+1)]);
    
    ind = find(ratioV(2:end) == min(ratioV(2:end)));
    channels{ind+1} = 'stim';

    inds = setdiff(inds,[1 (ind+1)]);

    ind = find(edges(inds,1) == min(edges(inds,1)));

    channels{inds(ind)} = 'air_puff';

    inds = setdiff(inds,inds(ind));

    ind = min(inds);
    channels{ind} = 'ch_a';
    ind = setdiff(inds,ind);
    channels{ind} = 'ch_b';
end


