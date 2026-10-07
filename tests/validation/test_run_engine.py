"""run_engine keeps every metabolite under its own name when the backend lacks one of them."""

import numpy as np

from basisremy.core.basisremy import BasisREMY
from tests.validation.cross_engine import run_engine


def test_run_engine_keeps_names_when_one_is_missing(monkeypatch):
    # a metabolite the backend lacks ('Mac' in FSL-MRS) must not shift the names after it
    br = BasisREMY('FSL-MRS')
    seen = {}

    def fake(params, *a, **k):
        seen['keys'] = list(params['Metabolites'])
        return {m: np.full(4, i + 1, dtype=complex) for i, m in enumerate(params['Metabolites'])}

    monkeypatch.setattr(br.backends['FSL-MRS'], 'run_simulation', fake)
    out, _ = run_engine(br, 'FSL-MRS', 'PRESS', 30.0, ['Lac', 'Mac', 'NAA', 'NAAG'], 3.0, 4, 2000.0)
    assert seen['keys'] == ['Lac', 'NAA', 'NAAG']
    assert {m: int(f[0].real) for m, f in out.items()} == {'Lac': 1, 'NAA': 2, 'NAAG': 3}
