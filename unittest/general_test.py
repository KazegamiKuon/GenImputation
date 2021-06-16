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

    def test_something(self):
        with g.reading('/home/cuong/VBDI/HungProject/GenImputation/data/raw/infiniumomni2-5-8-v1-3-a2.csv.gz') as ohg:
            for line in ohg:
                if line.startswith('[Assay]'):
                    break
            items = g.line_to_items(ohg.readline(),split=',')
            header_dict = g.list_to_dict(items)
            last_items = items
            for line in ohg:
                items = g.line_to_items(line,split=',')
                last_items = items
                pass
                # genome_build = items[header_dict['GenomeBuild']]
                # if genome_build == '38':
                #     print(genome_build)
                #     break        

if __name__ == '__main__':
    unittest.main()