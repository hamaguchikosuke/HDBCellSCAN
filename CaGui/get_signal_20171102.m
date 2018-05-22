function get_signal_20171102(ops1)

neuropilSub    = 'surround'; %getOr(ops1, {'neuropilSub'}, 'surround')

switch neuropilSub   
    case 'surround'
        get_signals_and_neuropil_kh(ops1, ops1.PlaneChString);
    otherwise
        error('Unknown option!');
end

