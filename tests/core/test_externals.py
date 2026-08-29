####################################################################################################
#                                         test_externals.py                                        #
####################################################################################################
#                                                                                                  #
# Purpose: Tests for core/externals.py's Octave patch step: a patch is applied once, is a no-op    #
#          on a patched tree, completes a tree that carries an older version of it file by file    #
#          (cached CI checkouts, upgrading users), and refuses a tree it does not fit.             #
#                                                                                                  #
####################################################################################################

from __future__ import annotations

import os
import subprocess
import sys

import pytest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

from basisremy.core import externals


def _git(repo, *args):
    return subprocess.run(['git', '-C', str(repo), *args], check=True, capture_output=True,
                          text=True, env={**os.environ, 'GIT_AUTHOR_NAME': 't', 'GIT_AUTHOR_EMAIL': 't@t',
                                          'GIT_COMMITTER_NAME': 't', 'GIT_COMMITTER_EMAIL': 't@t'}).stdout


@pytest.fixture
def repo(tmp_path):
    """A tiny git repo (two files) plus a two-file patch against it."""
    r = tmp_path / 'ext'
    r.mkdir()
    _git(r, 'init', '-q')
    (r / 'a.m').write_text('x=any(v,\'all\');\n')
    (r / 'b.m').write_text('y=sum(w,\'all\');\n')
    _git(r, 'add', '.')
    _git(r, 'commit', '-q', '-m', 'pinned')
    (r / 'a.m').write_text('x=any(any(v));\n')
    (r / 'b.m').write_text('y=sum(sum(w));\n')
    patch = tmp_path / 'octave.patch'
    patch.write_text(_git(r, 'diff'))
    _git(r, 'checkout', '--', '.')          # back to the pristine pinned tree
    return r, patch


def _with_patch(monkeypatch, patch):
    monkeypatch.setitem(externals.PATCHES, 'toy', str(patch))
    monkeypatch.setitem(externals.REGISTRY, 'toy', ('https://example.invalid/toy.git', 'deadbeefcafe'))


def test_apply_then_noop(repo, monkeypatch):
    r, patch = repo
    _with_patch(monkeypatch, patch)
    externals._apply_patch('toy', r, os.environ)
    assert (r / 'a.m').read_text() == 'x=any(any(v));\n'
    assert (r / 'b.m').read_text() == 'y=sum(sum(w));\n'
    externals._apply_patch('toy', r, os.environ)   # second call: nothing to do, no error
    assert _git(r, 'status', '--short').count('M') == 2


def test_completes_a_partially_patched_tree(repo, monkeypatch):
    r, patch = repo
    _with_patch(monkeypatch, patch)
    (r / 'a.m').write_text('x=any(any(v));\n')     # an older patch had touched a.m only
    externals._apply_patch('toy', r, os.environ)
    assert (r / 'b.m').read_text() == 'y=sum(sum(w));\n'
    assert (r / 'a.m').read_text() == 'x=any(any(v));\n'


def test_refuses_a_foreign_tree(repo, monkeypatch):
    r, patch = repo
    _with_patch(monkeypatch, patch)
    (r / 'b.m').write_text('y=something_else;\n')
    with pytest.raises(externals.ExternalFetchError, match='b.m'):
        externals._apply_patch('toy', r, os.environ)


def test_split_patch(repo):
    _, patch = repo
    chunks = externals._split_patch(patch.read_text())
    assert len(chunks) == 2 and chunks[0].startswith('diff --git a/a.m') and chunks[1].startswith('diff --git a/b.m')
