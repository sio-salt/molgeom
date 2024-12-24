import numpy as np
from numpy.typing import ArrayLike
from numpy.testing import assert_allclose


def xyz_to_pol(
    x: ArrayLike, y: ArrayLike, z: ArrayLike
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Calculates origin-centered spherical coordinates from Cartesian coordinates.

    Args:
        x, y, z: ArrayLike (e.g., lists, NumPy arrays)
    Returns:
        r, theta, phi: NumPy arrays
    """
    x, y, z = np.asarray(x), np.asarray(y), np.asarray(z)

    # Tolerance for error when checking near-zero values
    TH_ERROR = 10**-8

    # Calculate radius
    r = np.asarray(np.sqrt(x**2 + y**2 + z**2))

    # Avoid division by zero for theta
    with np.errstate(invalid="ignore", divide="ignore"):
        theta = np.asarray(
            np.arccos(np.clip(z / r, -1.0, 1.0))
        )  # Clip to avoid numerical issues

    # Calculate phi using arctan2
    phi = np.asarray(np.arctan2(y, x))
    if np.isscalar(phi):
        if phi < 0:
            phi += 2 * np.pi
    else:
        phi[phi < 0] += 2 * np.pi  # Ensure phi is in the range [0, 2π]

    # Handle special cases
    is_zero = (abs(x) < TH_ERROR) & (abs(y) < TH_ERROR) & (abs(z) < TH_ERROR)
    phi[is_zero] = -np.pi  # Undefined phi
    theta[is_zero] = -np.pi  # Undefined theta

    # Convert theta and phi to degrees
    theta = np.degrees(theta)
    phi = np.degrees(phi)

    return r, theta, phi


def pol_to_xyz(
    r: ArrayLike, theta: ArrayLike, phi: ArrayLike
) -> tuple[ArrayLike, ArrayLike, ArrayLike]:
    """
    args:
        r, theta, phi: ArrayLike
    return:
        x, y, z: ArrayLike
    """
    r, theta, phi = np.asarray(r), np.asarray(theta), np.asarray(phi)
    x = r * np.sin(np.radians(theta)) * np.cos(np.radians(phi))
    y = r * np.sin(np.radians(theta)) * np.sin(np.radians(phi))
    z = r * np.cos(np.radians(theta))
    return x, y, z


def test_xyz_to_pol():
    x = np.array([1, 0, 0, 0.628973895451, 0.619418368683, 0])
    y = np.array([0, 1, -0.628973895451, 0, 0.109220170745, 0])
    z = np.array([0, 0, -1.292593551292, -1.292593551292, -1.292593551292, 1])

    r, theta, phi = xyz_to_pol(x, y, z)
    assert_allclose(r, [1, 1, 1.4375, 1.4375, 1.4375, 1], atol=1e-6)
    assert_allclose(
        theta,
        [90, 90, 154.05252179901214, 154.05252179901214, 154.052521798993, 0],
        atol=1e-6,
    )
    assert_allclose(phi, [0, 90, 270, 0, 9.999999999978812, 0], atol=1e-6)

    del x, y, z, r, theta, phi

    x, y, z = 0, 0, 0
    r, theta, phi = xyz_to_pol(x, y, z)
    assert_allclose(r, 0, atol=1e-6)
    assert_allclose(theta, -180, atol=1e-6)
    assert_allclose(phi, -180, atol=1e-6)


if __name__ == "__main__":
    test_xyz_to_pol()
    print("All tests pass")
