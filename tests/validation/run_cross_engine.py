"""Run the cross-engine comparison and print markdown tables.

    BASISREMY_NO_DOCKER=1 python -m tests.validation.run_cross_engine press
    BASISREMY_NO_DOCKER=1 python -m tests.validation.run_cross_engine mega
    python -m tests.validation.run_cross_engine press FidaIdeal,FSL-MRS,Spant   # engine subset
"""

import sys

from basisremy.core.basisremy import BasisREMY
from tests.validation.cross_engine import PANEL, run_engine, agreement_table, format_table

B0, BW, N = 3.0, 2000.0, 2048


def main(which, engines=None):
    br = BasisREMY('FidaIdeal')
    results, timings = {}, {}
    if which == 'press':
        engines = engines or ['FidaIdeal', 'FSL-MRS', 'Spinach', 'Spant', 'Vespa']
        for e in engines:
            try:
                results[e], timings[e] = run_engine(br, e, 'PRESS', 35.0, PANEL, B0, N, BW)
            except Exception as exc:  # noqa: BLE001 - report and carry on
                print(f"!! {e} failed: {exc}")
        ref = engines[0]
    else:
        metabs = ['NAA', 'Cr', 'Glu', 'GABA', 'GSH']
        engines = engines or ['FidaMegaPressIdeal', 'FSL-MRS', 'Spant']
        for e in engines:
            try:
                results[e], timings[e] = run_engine(br, e, 'MEGA-PRESS', 68.0, metabs, B0, N, BW)
            except Exception as exc:  # noqa: BLE001
                print(f"!! {e} failed: {exc}")
        ref = engines[0]
    print(f"\n# {which}: reference {ref}\n")
    print(format_table(agreement_table(results, ref, BW, B0), timings))
    return results


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'press',
         sys.argv[2].split(',') if len(sys.argv) > 2 else None)
