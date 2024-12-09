from molgeom import FancyIndexingList


def test_FancyIndexing():
    fil = FancyIndexingList([1, 2, 3, 4, 5])

    assert fil[1] == 2
    assert fil[-1] == 5
    assert fil[[1, 2]] == [2, 3]
    assert fil[(1, 2)] == [2, 3]
    assert fil[1, 2] == [2, 3]
    assert fil[1, -2] == [2, 4]
    assert fil[[1, 2, 3]] == [2, 3, 4]
    assert fil[[0, -1, 2, -3]] == [1, 5, 3, 3]
    assert fil[1:] == [2, 3, 4, 5]
    assert fil[:3] == [1, 2, 3]
    assert fil[-2:] == [4, 5]
    assert fil[:-2] == [1, 2, 3]
    assert fil[1:3] == [2, 3]
    assert fil[:] == [1, 2, 3, 4, 5]
    assert fil[1:3:2] == [2]
    assert fil[1:4:2] == [2, 4]
    assert fil[1:5:2] == [2, 4]
    assert fil[1::2] == [2, 4]
    assert fil[:3:2] == [1, 3]
    assert fil[::2] == [1, 3, 5]
    assert fil[::-1] == [5, 4, 3, 2, 1]
    fil[1] = 10
    assert fil[1] == 10
    fil[[1, 2]] = [20, 30]
    assert fil[[1, 2]] == [20, 30]
    assert fil == [1, 20, 30, 4, 5]
    fil[1, 2, 3] = [20, 30, 40]
    assert fil[1, 2, 3] == [20, 30, 40]
    assert fil == [1, 20, 30, 40, 5]
    fil = FancyIndexingList([1, 2, 3, 4, 5])
    fil[[0, 1, 2, -1]] = [10, 20, 30, 40]
    assert fil[[0, 1, 2, -1]] == [10, 20, 30, 40]
    assert fil == [10, 20, 30, 4, 40]


def test_slice_id_check():
    fil = FancyIndexingList([1, 2, 3, 4, 5])

    assert fil[:] == fil and fil[:] is not fil
    assert fil[1:3] == FancyIndexingList([2, 3]) and fil[1:3] is not FancyIndexingList(
        [2, 3]
    )
