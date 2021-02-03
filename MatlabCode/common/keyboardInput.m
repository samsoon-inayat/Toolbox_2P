function occ = keyboardInput(cc,ccLims,steps,msg)
occ = cc;
if ~isempty(msg)
    display(msg);
end
keyVal = getkey;%pause;
if keyVal == 29 % Right arrow
    occ = cc + steps(1);
end
if keyVal == 28 % left arrow
    occ = cc - steps(1);
end
if keyVal == 30 % Up arrow
    occ = cc + steps(2);
end
if keyVal == 31 % down arrow
    occ = cc - steps(2);
end

if keyVal == 13 % l
    occ = input('Enter cell index : ');
end

if occ < ccLims(1)
    occ = ccLims(1);
end

if occ > ccLims(2)
    occ = ccLims(2);
end

if keyVal == 27 % Esc
    occ = -1;
end

