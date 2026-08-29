####################################################################################################
#                                     metabolite_identity.py                                       #
####################################################################################################
#                                                                                                  #
# Purpose: Say when two backend-native metabolite names mean the same compound, so a selection     #
#          made under one backend carries over to the next one. Nothing is renamed: engines and    #
#          exports keep their native keys; this table is only consulted when a selection moves     #
#          between backends (and can later back a uniform display naming).                         #
#                                                                                                  #
#          Rules:                                                                                  #
#            * names that differ only in case are one compound (spant writes 'naa', 'pcr', …);     #
#            * the explicit groups below join the remaining synonyms ('Ins' ≡ 'mI' ≡ 'ins');       #
#            * variant spin systems are NOT merged with their parent: 'GABA_gov', 'GSH_v2',        #
#              'Tau_govind', spant's '*_rt' / 'a_glc' / 'b_glc' / 'naa2', and FID-A's lumped       #
#              'Lip' versus spant's 'lip09'… stay separate entries.                                #
#                                                                                                  #
####################################################################################################

from __future__ import annotations

from collections.abc import Iterable

# Synonyms across the shipped backends (FID-A / FSL-MRS / Vespa / CustomSLaser share the
# first spelling, MRSCloud the second where present, spant the lower-case one).
IDENTITY_GROUPS = (
    ("Ins", "mI", "ins"),        # myo-inositol
    ("Scyllo", "sI", "sins"),    # scyllo-inositol
    ("PE", "peth"),              # phosphorylethanolamine
    ("bHG", "2hg"),              # 2-hydroxyglutarate
    ("Ch", "cho"),               # choline
)

_IDENTITY = {name.lower(): group[0].lower() for group in IDENTITY_GROUPS for name in group}


def metabolite_identity(name: str) -> str:
    """Backend-independent identity of a native metabolite name (lower-case key)."""
    key = str(name).lower()
    return _IDENTITY.get(key, key)


def same_metabolite(a: str, b: str) -> bool:
    return metabolite_identity(a) == metabolite_identity(b)


def translate_metabolites(names: Iterable[str], target_keys: Iterable[str]) -> list[str]:
    """Map native names of one backend onto the native keys of another.

    Names without a counterpart in `target_keys` are dropped; the order of
    `names` is kept.
    """
    by_identity = {metabolite_identity(k): k for k in target_keys}
    out = []
    for name in names:
        key = by_identity.get(metabolite_identity(name))
        if key is not None and key not in out:
            out.append(key)
    return out
