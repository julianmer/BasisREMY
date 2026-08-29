####################################################################################################
#                                   test_metabolite_identity.py                                    #
####################################################################################################
#                                                                                                  #
# Purpose: Tests for core/metabolite_identity.py — synonyms across backends are one compound,      #
#          variant spin systems stay separate, selections translate between backends, and no       #
#          backend carries two native keys with the same identity.                                 #
#                                                                                                  #
####################################################################################################

from __future__ import annotations

import importlib
import inspect
import os
import pkgutil
import sys

import pytest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

import basisremy.backends as backends_pkg
from basisremy.backends.base import Backend
from basisremy.core.metabolite_identity import (
    IDENTITY_GROUPS,
    metabolite_identity,
    same_metabolite,
    translate_metabolites,
)


def _all_backend_classes():
    out = []
    for m in pkgutil.iter_modules(backends_pkg.__path__):
        if m.name == 'base':
            continue
        mod = importlib.import_module(f'basisremy.backends.{m.name}')
        for _, cls in inspect.getmembers(mod, inspect.isclass):
            if cls.__module__ == mod.__name__ and issubclass(cls, Backend):
                out.append(cls)
    return out


@pytest.mark.parametrize("group", IDENTITY_GROUPS)
def test_groups_share_one_identity(group):
    idents = {metabolite_identity(n) for n in group}
    assert len(idents) == 1


@pytest.mark.parametrize("a,b", [
    ("Ins", "mI"), ("mI", "ins"), ("Scyllo", "sI"), ("sI", "sins"),
    ("PE", "peth"), ("bHG", "2hg"), ("Ch", "cho"),
    ("NAA", "naa"), ("PCr", "pcr"), ("GABA", "gaba"),   # case only
])
def test_synonyms(a, b):
    assert same_metabolite(a, b)


@pytest.mark.parametrize("a,b", [
    ("GABA", "GABA_gov"), ("GABA_gov", "GABA_govind"), ("GSH", "GSH_v2"),
    ("Tau", "Tau_govind"), ("gaba", "gaba_rt"), ("naa", "naa2"),
    ("glc", "a_glc"), ("Lip", "lip09"), ("Cr", "cr_ch3_rt"), ("PCh", "GPC"),
])
def test_variants_stay_separate(a, b):
    assert not same_metabolite(a, b)


def test_translate_fida_to_mrscloud_and_spant():
    fida = ['NAA', 'Ins', 'Scyllo', 'PE', 'Lip']
    assert translate_metabolites(fida, ['Asc', 'mI', 'NAA', 'sI', 'PE']) == ['NAA', 'mI', 'sI', 'PE']
    spant = ['naa', 'ins', 'sins', 'peth', '2hg', 'gaba_rt']
    assert translate_metabolites(fida, spant) == ['naa', 'ins', 'sins', 'peth']
    assert translate_metabolites(['bHG', 'GABA_gov'], spant) == ['2hg']


def test_translate_keeps_order_and_dedupes():
    assert translate_metabolites(['sI', 'NAA', 'Scyllo'], ['Scyllo', 'NAA']) == ['Scyllo', 'NAA']


@pytest.mark.parametrize("cls", _all_backend_classes(), ids=lambda c: c.__name__)
def test_backend_keys_have_unique_identities(cls):
    keys = list(cls().metabs)
    idents = [metabolite_identity(k) for k in keys]
    dupes = {i for i in idents if idents.count(i) > 1}
    assert not dupes, f"{cls.__name__}: several native keys share an identity: {dupes}"


def test_update_from_backend_carries_selection_across_names():
    from basisremy.backends.fida_backends import FidaIdeal
    from basisremy.backends.mrscloud_backend import MRSCloudBackend
    from basisremy.backends.spant_backend import SpantBackend

    fida = FidaIdeal()
    for k in fida.metabs:
        fida.metabs[k] = False
    fida.metabs['Ins'] = True
    fida.metabs['NAA'] = True
    fida.metabs['Scyllo'] = False

    mrs = MRSCloudBackend()
    mrs.metabs['sI'] = True          # must be switched off by the carried 'Scyllo' = False
    mrs.update_from_backend(fida)
    assert mrs.metabs['mI'] is True and mrs.metabs['NAA'] is True and mrs.metabs['sI'] is False
    assert 'mI' in mrs.mandatory_params['Metabolites']
    assert 'sI' not in mrs.mandatory_params['Metabolites']

    spant = SpantBackend()
    spant.update_from_backend(mrs)
    assert spant.metabs['ins'] is True and spant.metabs['naa'] is True and spant.metabs['sins'] is False
