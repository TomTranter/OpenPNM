import pytest
import numpy as np
import openpnm as op
import matplotlib.pyplot as plt
from numpy.testing import assert_allclose
from openpnm import topotools


class PlotToolsTest:

    def setup_class(self):
        self.ws = op.Workspace()

    def test_plot_tutorial(self):
        pn = op.network.Cubic(shape=[4, 4, 1])
        g = topotools.plot_tutorial(pn)
        plt.close()

    def test_plot_networkx(self):
        # 2D networks in XY, YZ, XZ planes
        for i in range(3):
            shape = np.ones(3, dtype=int)
            shape[np.arange(3) != i] = [5, 8]
            pn = op.network.Cubic(shape=shape)
            x, y = pn["pore.coords"].T[op.topotools.dimensionality(pn)]
            fig, ax = plt.subplots()
            m = op.topotools.plot_networkx(pn, ax=ax)
            x_plot, y_plot = np.array(m.get_offsets()).T
            np.testing.assert_allclose(x_plot, x)
            np.testing.assert_allclose(y_plot, y)
            plt.close()
        # 1D networks in XY, YZ, XZ planes
        for i in range(3):
            shape = np.ones(3, dtype=int)
            shape[np.arange(3) == i] = [5]
            pn = op.network.Cubic(shape=shape)
            x, = pn["pore.coords"].T[op.topotools.dimensionality(pn)]
            fig, ax = plt.subplots()
            m = op.topotools.plot_networkx(pn, ax=ax)
            x_plot, y_plot = np.array(m.get_offsets()).T
            np.testing.assert_allclose(x_plot, x)
            plt.close()

    def test_plot_networkx_3d(self):
        pn = op.network.Cubic(shape=[5, 8, 3])
        with pytest.raises(Exception):
            op.topotools.plot_networkx(pn)


if __name__ == '__main__':

    t = PlotToolsTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: ' + item)
            t.__getattribute__(item)()
