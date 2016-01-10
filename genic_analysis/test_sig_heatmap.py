import unittest

from sig_heatmap import meth_dict, stats_points, scores

class TestSigHeatmap(unittest.TestCase):
    def setUp(self):
        self.meth_data = "/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/tissue/Tissue/CG_all_4.txt"
        self.m = meth_dict(self.meth_data)

    def test_stats_points(self):
        mmax, mtypes, mtype = stats_points(self.m)
        self.assertEqual(mmax, 300)
        pos = [p for p,v in self.m[mtype]]
        pos.sort()
        ### upstream reagion is grouped by 20 and contains 100 regions
        self.assertEqual((pos[0],pos[99]),(-2000,-20))
        ### genebody grouped by 10bps and contains 100 regions
        self.assertEqual((pos[100],pos[199]),(0,990))
        ### downstream grouped by 20bps and contains 100 regions...
        ### numbering should not start at 1100 when making graphs need to change to 1000
        self.assertEqual((pos[200],pos[299]),(1100,3080))

    def test_scores(self):
        mmax, mtypes, mtype = stats_points(self.m)
        x =  range(0,mmax)
        ### capture gene body
        self.assertEqual(sorted(self.m[mtype])[100][0], 0)
        self.assertEqual(sorted(self.m[mtype])[199][0], 990)
        ##scores(mmax,20,mtypes,mtype,self.m)



if __name__ == '__main__':
    unittest.main()
    suite = unittest.TestLoader().loadTestsFromTestCase(TestSequenceFunctions)
    unittest.TextTestRunner(verbosity=2).run(suite)
