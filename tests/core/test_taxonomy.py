####################################################################################################
#                                       test_taxonomy.py                                           #
####################################################################################################
#                                                                                                  #
# Authors: J. P. Merkofer (j.p.merkofer@tue.nl)                                                    #
#                                                                                                  #
# Created: 25/04/26                                                                                #
#                                                                                                  #
# Purpose: Verifies the Category → Backend taxonomy that drives the new GUI cascade.               #
#                                                                                                  #
####################################################################################################

from __future__ import annotations

import os
import sys

import pytest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

from core.basisremy import BasisREMY


@pytest.fixture(scope='module')
def br():
    return BasisREMY()


def test_categories_match_expected(br):
    expected = {'FID-A', 'Custom', 'MRSCloud', 'FSL-MRS'}
    actual = {c for c, names in br.categories.items() if names}
    assert expected.issubset(actual), \
        f"Missing categories: {expected - actual}"


def test_every_backend_has_category(br):
    for name, inst in br.backends.items():
        assert inst.category, f"Backend {name!r} missing category"
        assert inst.category in br.categories, \
            f"Backend {name!r} category {inst.category!r} not registered"


def test_lcmodel_lives_under_fida(br):
    assert 'FidaIdeal' in br.categories['FID-A']
    assert br.backends['FidaIdeal'].category == 'FID-A'


def test_slaser_lives_under_custom(br):
    assert 'CustomSLaser' in br.categories['Custom']
    assert br.backends['CustomSLaser'].category == 'Custom'


def test_fida_category_has_multiple_options(br):
    """FID-A is the showcase category — must list >1 backend (Ideal + shaped)."""
    assert len(br.categories['FID-A']) > 1, \
        "FID-A category should contain Ideal + shaped sims"


def test_set_category_switches_backend(br):
    br.set_category('Custom')
    assert br.backend.category == 'Custom'
    br.set_category('FID-A')
    assert br.backend.category == 'FID-A'


def test_legacy_lcmodel_alias(br):
    """Old 'LCModel' identifier is gone — must raise ValueError."""
    with pytest.raises(ValueError):
        br.set_backend('LCModel')


def test_set_default_backend(br):
    br.set_backend('FidaIdeal')
    assert br.backend.name == 'FidaIdeal'


def test_unknown_category_raises(br):
    with pytest.raises(ValueError):
        br.set_category('NopeSim')


def test_display_names_present(br):
    """Every backend must expose a non-empty display_name for the GUI."""
    for inst in br.backends.values():
        assert inst.display_name, f"{inst.name} missing display_name"



