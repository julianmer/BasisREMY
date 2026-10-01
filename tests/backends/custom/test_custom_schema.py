"""CustomSLaser parseREMY: header voxel sizes (mm) reach the simulation in cm."""

from basisremy.backends.custom_backends import CustomSLaser


def test_parse_remy_voxel_mm_to_cm():
    b = CustomSLaser()
    mandatory, _ = b.parseREMY({'LeftRightSize': 30, 'AnteriorPosteriorSize': 20.0,
                                'ExcitationFlipAngle': 90})
    assert mandatory['thkX'] == 3.0
    assert mandatory['thkY'] == 2.0
    # FOV and the refocusing flip angle are not header values: None keeps the defaults
    assert mandatory['fovX'] is None and mandatory['fovY'] is None
    assert mandatory['Flip Angle'] is None


def test_parse_remy_missing_sizes_stay_blank():
    mandatory, _ = CustomSLaser().parseREMY({})
    assert mandatory['thkX'] is None and mandatory['thkY'] is None
