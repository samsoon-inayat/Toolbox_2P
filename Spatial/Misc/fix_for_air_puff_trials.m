function tei = fix_for_air_puff_trials(tei,dcfilename)

if ~exist('dcfilename','var')
    dcfilename = 'define_contexts.m';
end

cvs = get_define_contexts_file(tei,dcfilename);
if isfield(cvs,'air_puff_trials')
    tei.b.air_puff_r = tei.b.air_puff_r(cvs.air_puff_trials);
    tei.b.air_puff_f = tei.b.air_puff_f(cvs.air_puff_trials);
end

if isfield(cvs,'air_puff_f_trials')
    tei.b.air_puff_f = tei.b.air_puff_f(cvs.air_puff_f_trials);
end

if isfield(cvs,'air_puff_r_trials')
    tei.b.air_puff_r = tei.b.air_puff_r(cvs.air_puff_r_trials);
end

n = 0;