import sys
import pandas as pd
from lib.utils import general as g
import os
import unittest
import numpy as np

class GeneralTest(unittest.TestCase):

    def test_concat_dataframe(self):
        group_data = g.get_group_file('/home/cuong/VBDI/HungProject/GenImputation/data/test/paper_hybrid_result',key=lambda x: x.split('_')[:-2])
        paths = []
        for k, data in group_data[0]:
            paths.extend(list(data))
        data = g.concat_dataframe(paths,header=None,sep=' ')
        pass_nb_row = 0
        for i, path in enumerate(paths):
            temp = pd.read_csv(path,header=None,sep=' ')
            nb_row = temp.values.shape[0]
            x = np.all([data.values[pass_nb_row:pass_nb_row+nb_row] == temp.values])
            assert x == True, 'not concat'
            pass_nb_row += nb_row

if __name__ == '__main__':
    unittest.main()