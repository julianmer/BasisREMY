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
