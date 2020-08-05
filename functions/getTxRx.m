function [tx, rx] = getTxRx(bss, uts, bsIdx, utIdx, dataDirection)

switch(dataDirection)
    case 'DL'
        tx = bss([bss.nodeIdx] == bsIdx);
        rx = uts([uts.nodeIdx] == utIdx);
    case 'UL'
        tx = uts([uts.nodeIdx] == bsIdx);
        rx = bss([bss.nodeIdx] == utIdx);
    otherwise
        error("Data direction '%s' not recognized",params.dataDirection)
end

end