function md = get_mdataxml(filename)
xDoc = xmlread(filename);
    tt = xDoc.getElementsByTagName('stimuli');
    tt1 = tt.item(0);
    tt2 = tt1.getElementsByTagName('stimulus');
    for ss = 1:tt2.getLength
        tt3 = tt2.item(ss-1);
        temp = char(tt3.getAttribute('name'));
        db(ii).stimD(ss).name = temp;
        temp = char(tt3.getAttribute('units'));
        db(ii).stimD(ss).units = temp;
        if strcmp(temp,'secs')
            temp = char(tt3.getAttribute('duration'));
            db(ii).stimD(ss).duration = temp;
        end
        if strcmp(temp,'pulses')
            temp = char(tt3.getAttribute('distance'));
            db(ii).stimD(ss).distance = temp;
        end
    end

    tt = xDoc.getElementsByTagName('channels');
    tt1 = tt.item(0);
    tt2 = tt1.getElementsByTagName('channel');
    for ss = 1:tt2.getLength
        tt3 = tt2.item(ss-1);
        tempi = str2num(tt3.getAttribute('value'));
        temp = char(tt3.getAttribute('name'));
        md.channel{tempi} = temp;
    end