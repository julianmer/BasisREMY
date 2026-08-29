####################################################################################################
#                                          test_parallel.py                                        #
####################################################################################################
#                                                                                                  #
# Purpose: Tests for Backend.simulate_in_parallel — metabolites simulated concurrently over        #
#          several Octave processes: order of the result, progress counting, stop events, error    #
#          propagation, the worker count, and (with Docker) that DockerOctave.feval is reentrant   #
#          and that a parallel FID-A / Spinach run equals the sequential one.                      #
#                                                                                                  #
####################################################################################################

from __future__ import annotations

import os
import sys
import threading
import time

import numpy as np
import pytest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

from basisremy.backends.base import Backend


class _Fake(Backend):
    """A backend whose 'Octave' is a plain object; no runtime needed."""

    def __init__(self):
        super().__init__()
        self.octave = object()


def _docker_available():
    try:
        import docker
        docker.from_env().ping()
        return True
    except Exception:
        return os.path.exists(os.path.expanduser('~/.orbstack/run/docker.sock'))


class TestSimulateInParallel:

    def test_results_keep_item_order(self, monkeypatch):
        monkeypatch.setenv('BASISREMY_OCTAVE_WORKERS', '3')
        b = _Fake()
        seen = []

        def work(session, item):
            time.sleep(0.02 * (5 - item))      # later items finish first
            seen.append(item)
            return item * 10

        out = b.simulate_in_parallel([1, 2, 3, 4], work)
        assert list(out) == [1, 2, 3, 4] and out == {1: 10, 2: 20, 3: 30, 4: 40}
        assert sorted(seen) == [1, 2, 3, 4]

    def test_progress_counts_completions(self, monkeypatch):
        monkeypatch.setenv('BASISREMY_OCTAVE_WORKERS', '2')
        calls = []
        _Fake().simulate_in_parallel(['a', 'b', 'c'], lambda s, m: m,
                                     progress_callback=lambda i, n: calls.append((i, n)))
        assert calls == [(1, 3), (2, 3), (3, 3)]

    def test_stop_event_stops_submitting(self, monkeypatch):
        monkeypatch.setenv('BASISREMY_OCTAVE_WORKERS', '1')
        stop = threading.Event()

        def work(session, item):
            if item == 'b':
                stop.set()
            return item

        out = _Fake().simulate_in_parallel(['a', 'b', 'c', 'd'], work, stop_event=stop)
        assert 'a' in out and 'd' not in out

    def test_first_error_propagates(self, monkeypatch):
        monkeypatch.setenv('BASISREMY_OCTAVE_WORKERS', '2')

        def work(session, item):
            if item == 2:
                raise RuntimeError('boom')
            return item

        with pytest.raises(RuntimeError, match='boom'):
            _Fake().simulate_in_parallel([1, 2, 3], work)

    def test_worker_count(self, monkeypatch):
        monkeypatch.setenv('BASISREMY_OCTAVE_WORKERS', '7')
        assert Backend.octave_workers() == 7
        monkeypatch.setenv('BASISREMY_OCTAVE_WORKERS', 'nonsense')
        assert 1 <= Backend.octave_workers() <= 4
        monkeypatch.delenv('BASISREMY_OCTAVE_WORKERS')
        assert 1 <= Backend.octave_workers() <= 4

    def test_sessions_are_shared_across_threads(self, monkeypatch):
        monkeypatch.setenv('BASISREMY_OCTAVE_WORKERS', '3')
        b = _Fake()
        sessions = set()
        _Fake.simulate_in_parallel(b, range(6), lambda s, m: sessions.add(id(s)) or m)
        assert sessions == {id(b.octave)}     # the (reentrant) runner is shared, not copied


@pytest.mark.skipif(not _docker_available(), reason='needs Docker')
class TestParallelLive:

    def test_docker_feval_is_reentrant(self):
        from basisremy.core.octave_manager import OctaveManager
        octave = OctaveManager().initialize_octave(prefer_docker=True)
        results = {}

        def call(k):
            results[k] = float(octave.feval('sum', [1.0, float(k)]))

        threads = [threading.Thread(target=call, args=(k,)) for k in range(4)]
        for th in threads:
            th.start()
        for th in threads:
            th.join()
        assert results == {0: 1.0, 1: 2.0, 2: 3.0, 3: 4.0}

    @pytest.mark.parametrize('backend', ['FidaIdeal', 'Spinach'])
    def test_parallel_equals_sequential(self, backend, monkeypatch):
        from basisremy.core.basisremy import BasisREMY
        params = {'Sequence': 'PRESS', 'Samples': 2048, 'Bandwidth': 2000, 'Bfield': 3.0,
                  'Linewidth': 1, 'TE': 35, 'TE2': 0, 'TM': 10, 'Center Freq': 127.7,
                  'Metabolites': ['NAA', 'Cr', 'Lac', 'Glu']}
        monkeypatch.setenv('BASISREMY_OCTAVE_WORKERS', '1')
        seq = BasisREMY().backends[backend].run_simulation(dict(params))
        monkeypatch.setenv('BASISREMY_OCTAVE_WORKERS', '4')
        par = BasisREMY().backends[backend].run_simulation(dict(params))
        assert list(par) == list(seq) == params['Metabolites']
        for m in seq:
            np.testing.assert_allclose(par[m], seq[m], rtol=0, atol=1e-12)
