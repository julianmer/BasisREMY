####################################################################################################
#                                       test_gui_smoke.py                                          #
####################################################################################################
#                                                                                                  #
# Authors: J. P. Merkofer (j.p.merkofer@tue.nl)                                                    #
#                                                                                                  #
# Created: 27/08/26                                                                                #
#                                                                                                  #
# Purpose: Browser-free smoke tests of the NiceGUI app via nicegui.testing.User: the page builds,  #
#          the three-step flow navigates, backends switch, and Simulate stays locked until the     #
#          visible parameters are filled.                                                          #
#                                                                                                  #
####################################################################################################


#*************#
#   imports   #
#*************#
import sys, os
import pytest
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))

from nicegui import ui
from nicegui.testing import User

from basisremy.gui.application import build_page


@pytest.fixture(autouse=True)
def _register_page(user: User):
    # must depend on `user` so registration happens after its globals reset
    ui.page('/')(build_page)
    yield


async def _open(user: User) -> None:
    """Open the page and wait past the start-up splash."""
    await user.open('/')
    await user.should_see('Skip', retries=100)


async def test_page_builds_with_dropzone(user: User) -> None:
    await _open(user)
    await user.should_see('Drop MRS data or click to browse')
    await user.should_see('Skip')


async def test_skip_reaches_parameters(user: User) -> None:
    await _open(user)
    user.find('Skip').click()
    await user.should_see('Engine')
    await user.should_see('Sequence file')
    await user.should_see('Metabolites')
    # pulses live in the sequence designer only
    await user.should_not_see('Excitation')
    await user.should_not_see('Refocusing')


async def test_simulate_locked_when_params_blank(user: User) -> None:
    await _open(user)
    user.find('Skip').click()
    await user.should_see('Simulate basis set')
    button = user.find('Simulate basis set').elements.pop()
    assert not button.enabled, \
        "Simulate must stay locked while scan parameters are blank"


async def test_continue_disabled_without_file(user: User) -> None:
    await _open(user)
    button = user.find('Continue').elements.pop()
    assert not button.enabled


async def test_metab_buttons_all_none_default(user: User) -> None:
    await _open(user)
    user.find('Skip').click()
    await user.should_see('Metabolites')
    user.find('None').click()
    await user.should_see('Metabolites')
    user.find('All').click()
    await user.should_see('Metabolites')
    user.find('Default').click()
    await user.should_see('Metabolites')


async def _select(user, pred, kind=ui.select, tries=50):
    import asyncio
    for _ in range(tries):
        try:
            found = [s for s in user.find(kind).elements if pred(s)]
        except AssertionError:          # none of that kind on the page yet
            found = []
        if found:
            return found[0]
        await asyncio.sleep(0.05)
    raise AssertionError(f'{kind.__name__} not found')


async def _engine(user, category):
    """Switch the engine and return the rebuilt sheet's sequence picker."""
    await user.should_see('Engine')
    engine = await _select(user, lambda s: 'FID-A' in s.options and 'MRSCloud' in s.options)
    engine.value = category
    await _select(user, lambda s: 'FID-A' in s.options and s.value == category and s is not engine)
    return await _select(user, lambda s: 'Spin Echo train' in s.options)


async def test_sequence_panel_opens(user: User) -> None:
    await _open(user)
    user.find('Skip').click()
    await user.should_see('Engine')
    user.find(marker='sequence-panel').click()
    await user.should_see('Sequence designer')


async def test_unsupported_sequences_are_greyed(user: User) -> None:
    await _open(user)
    user.find('Skip').click()
    seq = await _engine(user, 'Vespa')
    disabled = {seq._values[i] for i, o in enumerate(seq.props['options']) if o.get('disable')}
    assert 'sLASER' in disabled and 'MEGA-PRESS' in disabled
    assert 'PRESS' not in disabled and 'STEAM' not in disabled


async def test_designer_saves_and_selects_the_design(user: User, tmp_path, monkeypatch) -> None:
    # Vespa PRESS: a standard refocusing pulse in the designer, Save -> the .seq is selected in the
    # sheet's file field and Vespa runs it as 'PRESS shaped' with the pulse from the file
    monkeypatch.setenv('BASISREMY_SEQUENCES_DIR', str(tmp_path))
    await _open(user)
    user.find('Skip').click()
    seq = await _engine(user, 'Vespa')
    seq.value = 'PRESS'
    await user.should_see('Echo split: symmetric (TE/2 each) in Vespa.')   # rebuilt PRESS sheet
    te = user.find(marker='param:TE').elements.pop()
    te.value = '30'
    user.find(marker='sequence-panel').click()
    await user.should_see('Sequence designer')
    refoc = await _select(user, lambda s: 'standard:sinc-ref' in s.options and s.value == 'ideal')
    refoc.value = 'standard:sinc-ref'
    await _select(user, lambda s: 'standard:sinc-ref' in s.options
                  and s.value == 'standard:sinc-ref' and s is not refoc)
    await user.should_see('Duration [ms]')
    user.find(marker='designer-save').click()
    saved = await _select(user, lambda s: any(str(k).endswith('.seq') for k in s.options)
                          and str(s.value).endswith('.seq'))
    assert str(saved.value).startswith(str(tmp_path))
    await user.should_see('Runs as given.')
    # the timings now come from the file; MRSCloud runs a PRESS design too
    assert 'readonly' in user.find(marker='param:TE').elements.pop().props
    engine = await _select(user, lambda s: 'FID-A' in s.options and 'MRSCloud' in s.options)
    assert 'cannot run' not in engine.options['MRSCloud']


async def test_echo_split_recommended_from_te(user: User) -> None:
    # FID-A ideal PRESS: typing TE fills the echo split (TE/2) as a recommendation
    await _open(user)
    user.find('Skip').click()
    seq = await _engine(user, 'FID-A')
    seq.value = 'PRESS'
    await user.should_see('TE2 (second echo) [ms]')
    te = user.find(marker='param:TE').elements.pop()
    assert 'br-v-missing' in te.classes
    te.value = '30'
    split = user.find(marker='param:TE2').elements.pop()
    assert split.value == '15' and 'br-v-rec' in split.classes
    assert 'br-v-user' in te.classes
