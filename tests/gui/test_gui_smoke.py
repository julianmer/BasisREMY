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


async def test_page_builds_with_dropzone(user: User) -> None:
    await user.open('/')
    await user.should_see('Drop MRS data or click to browse')
    await user.should_see('Skip')


async def test_skip_reaches_parameters(user: User) -> None:
    await user.open('/')
    user.find('Skip').click()
    await user.should_see('Simulation Software')
    await user.should_see('Metabolites')


async def test_simulate_locked_when_params_blank(user: User) -> None:
    await user.open('/')
    user.find('Skip').click()
    await user.should_see('Simulate basis set')
    button = user.find('Simulate basis set').elements.pop()
    assert not button.enabled, \
        "Simulate must stay locked while scan parameters are blank"


async def test_continue_disabled_without_file(user: User) -> None:
    await user.open('/')
    button = user.find('Continue').elements.pop()
    assert not button.enabled


async def test_metab_buttons_all_none_default(user: User) -> None:
    await user.open('/')
    user.find('Skip').click()
    await user.should_see('Metabolites')
    user.find('None').click()
    await user.should_see('Metabolites')
    user.find('All').click()
    await user.should_see('Metabolites')
    user.find('Default').click()
    await user.should_see('Metabolites')


async def test_sequence_panel_opens(user: User) -> None:
    await user.open('/')
    user.find('Skip').click()
    await user.should_see('Sequence…')
    user.find('Sequence…').click()
    await user.should_see('Runs on')
    await user.should_see('Whole sequence from a file')
    await user.should_see('Pulseq (.seq)')


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


async def test_sequence_panel_switches_ideal_and_waveform(user: User) -> None:
    # Vespa PRESS: a standard refocusing pulse turns it into 'PRESS shaped', Ideal turns it back
    await user.open('/')
    user.find('Skip').click()
    await user.should_see('Simulation Software')
    software = next(s for s in user.find(ui.select).elements if s.value == 'MRSCloud')
    software.value = 'Vespa'
    seq = await _select(user, lambda s: 'PRESS shaped' in s.options)
    seq.value = 'PRESS'
    user.find('Sequence…').click()
    await user.should_see('Refocusing')
    refoc = await _select(user, lambda s: 'standard:sinc-ref' in s.options)
    assert refoc.value == 'ideal'
    refoc.value = 'standard:sinc-ref'
    seq = await _select(user, lambda s: 'PRESS shaped' in s.options and s.value == 'PRESS shaped')
    await _select(user, lambda n: n.props.get('label') == 'Duration [ms] (RefTp)', ui.number)
    refoc = await _select(user, lambda s: 'standard:sinc-ref' in s.options)
    assert refoc.value == 'standard:sinc-ref'
    refoc.value = 'ideal'
    seq = await _select(user, lambda s: 'PRESS shaped' in s.options and s.value == 'PRESS')
    with pytest.raises(AssertionError):     # no duration field once the pulse is ideal again
        await _select(user, lambda n: n.props.get('label') == 'Duration [ms] (RefTp)', ui.number, tries=5)
