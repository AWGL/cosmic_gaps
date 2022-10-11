import unittest
import pandas
from filter_raw_data import *


class test_filter_raw_data(unittest.TestCase):

    def test_prefilter_cosmic_data(self):


        '''
        check that the headings get capitalised and that the gene list is as expected
        '''


        with open("config.yaml") as f:
            config = yaml.load(f, Loader=yaml.FullLoader)
        file=prefilter_cosmic_data("test_data/test_cosmic_output.tsv", config)
        print(file)

        self.assertEqual(len(file),12)

	#check that the headers have all been capitalized and any spaces have been replaced with "_"

        gene_name=list(file["GENE_NAME"])
        accession_number=list(file["ACCESSION_NUMBER"])
        gene_cds_length=list(file["GENE_CDS_LENGTH"])
        hgnc_id=list(file["HGNC_ID"])
        sample_name=list(file["SAMPLE_NAME"])
        id_sample=list(file["ID_SAMPLE"])
        id_tumour=list(file["ID_TUMOUR"])
        primary_site=list(file["PRIMARY_SITE"])
        site_subtype_1=list(file["SITE_SUBTYPE_1"])
        site_subtype_2=list(file["SITE_SUBTYPE_2"])
        site_subtype3=list(file["SITE_SUBTYPE_3"])
        primary_histology=list(file["PRIMARY_HISTOLOGY"])
        histology_subtype_1=list(file["HISTOLOGY_SUBTYPE_1"])
        histology_subtype_2=list(file["HISTOLOGY_SUBTYPE_2"])
        histology_subtype_3=list(file["HISTOLOGY_SUBTYPE_3"])
        genomic_wide_screen=list(file["GENOME_WIDE_SCREEN"])
        genomic_mutation_id=list(file["GENOMIC_MUTATION_ID"])
        legacy_mutation_id=list(file["LEGACY_MUTATION_ID"])
        mutation_id=list(file["MUTATION_ID"])
        mutation_cds=list(file["MUTATION_CDS"])
        mutation_aa=list(file["MUTATION_AA"])
        mutation_description=list(file["MUTATION_DESCRIPTION"])
        mutation_zygosity=list(file["MUTATION_ZYGOSITY"])
        loh=list(file["LOH"])
        grch=list(file["GRCH"])
        mutation=list(file["MUTATION_GENOME_POSITION"])
        mutation_strand=list(file["MUTATION_STRAND"])
        snp=list(file["SNP"])
        resistance_mutation=list(file["RESISTANCE_MUTATION"])
        fathmm_prediction=list(file["FATHMM_PREDICTION"])
        fathmm_score=list(file["FATHMM_SCORE"])
        mutation_somatic_status=list(file["MUTATION_SOMATIC_STATUS"])
        pubmed_pmid=list(file["PUBMED_PMID"])
        id_study=list(file["ID_STUDY"])
        sample_type=list(file["SAMPLE_TYPE"])
        tumour_origin=list(file["TUMOUR_ORIGIN"])
        age=list(file["AGE"])
        hgvsp=list(file["HGVSP"])
        hgvsc=list(file["HGVSC"])
        hgvsg=list(file["HGVSG"])


	#check that the list is as expected- KRAS2 is in this list but that is not a problem if it gets removed from final output
        self.assertEqual(gene_name,["EGFR", "BRAF", "BRAF", "KRAS", "KRAS2", "NRAS", "EGFR", "EGFR", "BRAF", "KIT", "KIT", "PDGFRA"])



    def test_prefilter_cosmic_data(self):

        '''
        check that a file is produced for each of the referrral types in the config file
        '''


        with open("config.yaml") as f:
            config = yaml.load(f, Loader=yaml.FullLoader)


        df=prefilter_cosmic_data("test_data/test_cosmic_output.tsv", config)

        referral_file=filter_by_referral(df, config)

        
        #check it only filters the lung EGFR rows

        Lung_unsorted=pandas.read_csv("Lung.csv")

        Lung=Lung_unsorted.sort_values(by='CHR')
        chr=list(Lung["CHR"])
        start=list(Lung["START"])
        end=list(Lung["END"])
        gene_name=list(Lung["GENE_NAME"])
        hgvs_c=list(Lung["HGVS_C"])
        hgvs_p=list(Lung["HGVS_P"])
        hgvs_g=list(Lung["HGVS_G"])
        count=list(Lung["COUNT"])

        self.assertEqual(chr,[11,13])
        self.assertEqual(start,[5896748, 5896750])
        self.assertEqual(end,[5896749, 5896751])
        self.assertEqual(gene_name,["EGFR", "EGFR"])
        self.assertEqual(hgvs_c,["HGVSC", "HGVSC"])
        self.assertEqual(hgvs_p,["HGVSP", "HGVSP"])
        self.assertEqual(hgvs_g,["11:g.5896749C>T", "13:g.5896751C>T"])
        self.assertEqual(count,[1, 1])



        gist_unsorted=pandas.read_csv("GIST.csv")

        gist=gist_unsorted.sort_values(by='CHR')
        print(gist)

        chr=list(gist["CHR"])
        start=list(gist["START"])
        end=list(gist["END"])
        gene_name=list(gist["GENE_NAME"])
        hgvs_c=list(gist["HGVS_C"])
        hgvs_p=list(gist["HGVS_P"])
        hgvs_g=list(gist["HGVS_G"])
        count=list(gist["COUNT"])

        self.assertEqual(chr,[17,18])
        self.assertEqual(start,[5896754,5896755])
        self.assertEqual(end,[5896755,5896756])
        self.assertEqual(gene_name,["KIT","PDGFRA"])
        self.assertEqual(hgvs_c,["HGVSC", "HGVSC"])
        self.assertEqual(hgvs_p,["HGVSP", "HGVSP"])
        self.assertEqual(hgvs_g,["17:g.5896755C>T", "18:g.5896756C>T"])
        self.assertEqual(count,[1, 1])



        #check that KIT2 does not get added to this list
        skin_unsorted=pandas.read_csv("Skin.csv")

        skin=skin_unsorted.sort_values(by='CHR')
        print(skin)

        chr=list(skin["CHR"])
        start=list(skin["START"])
        end=list(skin["END"])
        gene_name=list(skin["GENE_NAME"])
        hgvs_c=list(skin["HGVS_C"])
        hgvs_p=list(skin["HGVS_P"])
        hgvs_g=list(skin["HGVS_G"])
        count=list(skin["COUNT"])

        self.assertEqual(chr,[14,15])
        self.assertEqual(start,[5896751,5896752])
        self.assertEqual(end,[5896752,5896753])
        self.assertEqual(gene_name,["BRAF","KIT"])
        self.assertEqual(hgvs_c,["HGVSC", "HGVSC"])
        self.assertEqual(hgvs_p,["HGVSP", "HGVSP"])
        self.assertEqual(hgvs_g,["14:g.5896752C>T", "15:g.5896753C>T"])
        self.assertEqual(count,[1, 1])


        #check the counts get added for overlapping regions

        colorectal_unsorted=pandas.read_csv("Colorectal.csv")

        colorectal=colorectal_unsorted.sort_values(by='CHR')
        chr=list(colorectal["CHR"])
        start=list(colorectal["START"])
        end=list(colorectal["END"])
        gene_name=list(colorectal["GENE_NAME"])
        hgvs_c=list(colorectal["HGVS_C"])
        hgvs_p=list(colorectal["HGVS_P"])
        hgvs_g=list(colorectal["HGVS_G"])
        count=list(colorectal["COUNT"])

        self.assertEqual(chr,[5,6,8])
        self.assertEqual(start,[5896741, 5896743, 5896745])
        self.assertEqual(end,[5896742, 5896744, 5896746])
        self.assertEqual(gene_name,["BRAF" ,"KRAS", "NRAS"])
        self.assertEqual(hgvs_c,["HGVSC", "HGVSC", "HGVSC"])
        self.assertEqual(hgvs_p,["HGVSP", "HGVSP", "HGVSP"])
        self.assertEqual(hgvs_g,["5:g.5896742C>T", "6:g.5896744C>T", "8:g.5896746C>T"])
        self.assertEqual(count,[2, 1, 1])




        #check for the breast referral type

        breast_unsorted=pandas.read_csv("Breast.csv")

        breast=breast_unsorted.sort_values(by='CHR')
        chr=list(breast["CHR"])
        start=list(breast["START"])
        end=list(breast["END"])
        gene_name=list(breast["GENE_NAME"])
        hgvs_c=list(breast["HGVS_C"])
        hgvs_p=list(breast["HGVS_P"])
        hgvs_g=list(breast["HGVS_G"])
        count=list(breast["COUNT"])

        self.assertEqual(chr,[18])
        self.assertEqual(start,[5896755])
        self.assertEqual(end,[5896756])
        self.assertEqual(gene_name,["PIK3CA"])
        self.assertEqual(hgvs_c,["HGVSC"])
        self.assertEqual(hgvs_p,["HGVSP"])
        self.assertEqual(hgvs_g,["18:g.5896756C>T"])
        self.assertEqual(count,[1])





