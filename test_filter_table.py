import unittest
from filter_table import *


class test_filter_table(unittest.TestCase):


    def test_filter_table(self):

        referral_list=['Melanoma', 'Lung', 'Colorectal', 'GIST']


        
        file=filter_table("Sample1", "Colorectal", "./test_data/", "./test_data/", referral_list )

        self.assertEqual(len(file),2)

        chr=list(file["Chr"])
        start=list(file["Start"])
        end=list(file["End"])
        gene=list(file["Gene"])
        info=list(file["Info"])
        counts=list(file["Counts"])
        percentage=list(file["Percentage"])

        self.assertEqual(chr[0], 1)
        self.assertEqual(start[0], 3)
        self.assertEqual(end[0], 15)
        self.assertEqual(gene[0], "gene1")
        self.assertEqual(info[0], "gene1(hgvs_c)")
        self.assertEqual(counts[0], 131)
        self.assertEqual(percentage[0], 13.83)

        self.assertEqual(chr[1], 1)
        self.assertEqual(start[1], 67)
        self.assertEqual(end[1], 87)
        self.assertEqual(gene[1], "gene3")
        self.assertEqual(info[1], "gene3(hgvs_c)")
        self.assertEqual(counts[1], 0)
        self.assertEqual(percentage[1], 0.0)


        file=filter_table("Sample2", "Melanoma", "./test_data/", "./test_data/", referral_list )

        self.assertEqual(len(file),2)

        chr=list(file["Chr"])
        start=list(file["Start"])
        end=list(file["End"])
        gene=list(file["Gene"])
        info=list(file["Info"])
        counts=list(file["Counts"])
        percentage=list(file["Percentage"])

        self.assertEqual(chr[0], 1)
        self.assertEqual(start[0], 3)
        self.assertEqual(end[0], 15)
        self.assertEqual(gene[0], "gene1")
        self.assertEqual(info[0], "gene1(hgvs_c)")
        self.assertEqual(counts[0], 133)
        self.assertEqual(percentage[0], 2660)


        self.assertEqual(chr[1], 1)
        self.assertEqual(start[1], 3)
        self.assertEqual(end[1], 20)
        self.assertEqual(gene[1], "gene1")
        self.assertEqual(info[1], "gene1(hgvs_c)")
        self.assertEqual(counts[1], 5)
        self.assertEqual(percentage[1], 100.0)



        file=filter_table("Sample3", "Lung", "./test_data/", "./test_data/", referral_list )

        self.assertEqual(len(file),5)


        chr=list(file["Chr"])
        start=list(file["Start"])
        end=list(file["End"])
        gene=list(file["Gene"])
        info=list(file["Info"])
        counts=list(file["Counts"])
        percentage=list(file["Percentage"])

        self.assertEqual(chr[0], 1)
        self.assertEqual(start[0], 3)
        self.assertEqual(end[0], 15)
        self.assertEqual(gene[0], "gene1")
        self.assertEqual(info[0], "gene1(hgvs_c)")
        self.assertEqual(counts[0], 0)
        self.assertEqual(percentage[0], 0.0)

        self.assertEqual(chr[1], 3)
        self.assertEqual(start[1], 7486)
        self.assertEqual(end[1], 8698)
        self.assertEqual(gene[1], "gene1")
        self.assertEqual(info[1], "gene1(hgvs_c)")
        self.assertEqual(counts[1], 0)
        self.assertEqual(percentage[1], 0.0)


        self.assertEqual(chr[2], 4)
        self.assertEqual(start[2], 947)
        self.assertEqual(end[2], 1957)
        self.assertEqual(gene[2], "gene1")
        self.assertEqual(info[2], "gene1(hgvs_c)")
        self.assertEqual(counts[2], 0)
        self.assertEqual(percentage[2],0.0)


        self.assertEqual(chr[3], 10)
        self.assertEqual(start[3], 22)
        self.assertEqual(end[3], 474)
        self.assertEqual(gene[3], "gene1")
        self.assertEqual(info[3], "gene1(hgvs_c)")
        self.assertEqual(counts[3], 0)
        self.assertEqual(percentage[3], 0.0)


        self.assertEqual(chr[4], 15)
        self.assertEqual(start[4], 285638)
        self.assertEqual(end[4], 385739)
        self.assertEqual(gene[4], "gene1")
        self.assertEqual(info[4], "gene1(hgvs_c)")
        self.assertEqual(counts[4], 0)
        self.assertEqual(percentage[4],0.0)



        file=filter_table("Sample4", "GIST", "./test_data/", "./test_data/", referral_list )

        self.assertEqual(len(file),3)

        chr=list(file["Chr"])
        start=list(file["Start"])
        end=list(file["End"])
        gene=list(file["Gene"])
        info=list(file["Info"])
        counts=list(file["Counts"])
        percentage=list(file["Percentage"])

        self.assertEqual(chr[0], 1)
        self.assertEqual(start[0], 3)
        self.assertEqual(end[0], 15)
        self.assertEqual(gene[0], "gene")
        self.assertEqual(info[0], "gene(hgvs_c)")
        self.assertEqual(counts[0], 0)
        self.assertEqual(percentage[0], 0.0)


        self.assertEqual(chr[1], 1)
        self.assertEqual(start[1], 3)
        self.assertEqual(end[1], 15)
        self.assertEqual(gene[1], "gene2")
        self.assertEqual(info[1], "gene2(hgvs_c)")
        self.assertEqual(counts[1], 19)
        self.assertEqual(percentage[1], 19.0)


        self.assertEqual(chr[2], 1)
        self.assertEqual(start[2], 999)
        self.assertEqual(end[2], 1011)
        self.assertEqual(gene[2], "gene2")
        self.assertEqual(info[2], "gene2(hgvs_c)")
        self.assertEqual(counts[2], 7)
        self.assertEqual(percentage[2],7)

