import pandas
import argparse
import decimal
from decimal import Decimal
import numpy

def calculate_percentage(x):

    '''
    function used in filter table function to calculate the percentage of cosmic variants in the gap compared to the whole gene

    '''

    if x['Counts'] == 0:
        percentage = 0
    else:
        percentage = (x['Counts']/x['Total_counts'])*100

    return percentage
    


def filter_table(sampleId, referral, gaps_path, bedfile_path, referral_list):

    '''
    uses the hotspots gaps file and intersect file to calculate the number of cosmic variants in each sequencing gap 

    '''

    #load file
    original_gaps_file = pandas.read_csv(gaps_path+sampleId+"_"+referral+"_hotspots.gaps", sep='\t', names=['Chr','Start', 'End', 'Info'])

    # only run loop when referral is in the referral list and the length of the gaps file is greater than 0
    if ((referral in referral_list) and (original_gaps_file.shape[0]!=0)):

        file = pandas.read_csv(gaps_path+sampleId+"_"+referral+"_intersect.txt", sep = '\t', names = ['Chr','Start', 'End', 'Info', 'Chr_cosmic','Start_cosmic', 'End_cosmic', '7', '8','9', '10', 'Counts', 'Overlap'])

        #For regions where there is no overlap make Counts value 0
        file['Counts'] = numpy.where(file['Overlap'] == 0, 0, file['Counts'])

        #only keep certain columns
        file[['Gene','Ignore']] = file.Info.str.split("(",expand=True,)
        file=file.filter(items = ['Chr', 'Start', 'End','Info', 'Gene', 'Counts'])

        #combine the rows for the same region and add the counts column for these rows
        file['Counts'] = file['Counts'].apply(lambda x: int(x))
        grouped_file = file.groupby(['Chr', 'Start', 'End', 'Info', 'Gene'], as_index = False).aggregate({'Counts': 'sum'})

        #Create a dictionary of the number of variants in the cosmic file for each of the genes 
        file = pandas.read_csv(bedfile_path+referral+".bed", sep='\t', names = ['Chr','Start', 'End', 'Gene', 'ENST','ENSP','Info', 'Counts'])
        genes = file['Gene']
        genes = list(set(genes))

        gene_counts_dict = {}
        for gene in genes:
            gene_counts=file[file["Gene"] == gene]
            gene_counts['Counts'] = gene_counts['Counts'].apply(lambda x: int(x))
            Total = gene_counts['Counts'].sum()
            gene_counts_dict[gene] = Total


        #Add the total counts and percentage column for each of the gaps
        grouped_file['Total_counts'] = grouped_file['Gene'].apply(lambda x: 0 if x not in genes else gene_counts_dict.get(x))
        grouped_file['Percentage'] = grouped_file.apply(lambda x: calculate_percentage(x), axis = 1)

        #Round the Percentage column to 2dp
        decimal.getcontext().rounding = decimal.ROUND_HALF_UP
        grouped_file['Percentage'] = grouped_file['Percentage'].apply(lambda x: float((Decimal(str(x)).quantize(Decimal('0.01')))))

        #output table to csv file
        grouped_file = grouped_file.filter(items = ['Chr', 'Start', 'End','Info', 'Gene', 'Counts', 'Percentage'])
        grouped_file.to_csv(gaps_path+sampleId+"_"+referral+"_output.csv", sep = ',', index = False)

    else:
        grouped_file = original_gaps_file
        grouped_file['Counts'] = 'N/A'
        grouped_file['Percentage'] = 'N/A'
        grouped_file.to_csv(gaps_path+sampleId+"_"+referral+"_output.csv", sep=',', index=False)

    return grouped_file



if __name__ == '__main__':


    parser=argparse.ArgumentParser()
    parser.add_argument('--sampleId', required=True)
    parser.add_argument('--referral', required=True)
    parser.add_argument('--gaps_path', required=True)
    parser.add_argument('--bedfile_path', required=True)
    args=parser.parse_args()

    sampleId=args.sampleId
    referral=args.referral
    gaps_path=args.gaps_path
    bedfile_path=args.bedfile_path
    referral_list=['Melanoma', 'Lung', 'Colorectal', 'GIST']


    if (referral!= "null"):
        filter_table(sampleId, referral, gaps_path, bedfile_path, referral_list)









