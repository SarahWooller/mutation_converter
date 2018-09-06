
import pandas as pd
import os
import numpy as np
import argparse


parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, epilog = """

The NetSurfP server predicts the surface accessibility and secondary structure of amino acids 
in an amino acid sequence. However, on occasion it is useful to find out how the surface accessibility and secondary structure changes when a gene
is mutated. 


In order to facilitate this, this function takes as its input a series of mutations stored as a text file, and outputs a series of files, with both the 
wild fasta and mutated fasta associated with the mutation ready for uploading to NetSurfP. 
Each file contains no more than 10,000 characters in order to not overwhelm the NetSurfP servers.


Usage:


Inputs:


path        -a path to text file of mutations of the form PIK3CA_E545K
gene_type   -either Gene_name or one of the other options listed above
output_dire -directory name where the results will be stored


Further information about the Inputs is given above.


Outputs:


Two files will be saved to the output directory of choice. One - wilds0 will give the wild fastas for each of the mutations.
Where several mutations share the same fastas the title of the fasta encorporates both mutations separated by a | as shown below.


>PIK3CA_E545K|PIK3CA_H1047R MPPRPSSGELWGIHLMPPRILVECLLPNGMIVTLECLREATLITIKHELFKEARKYPLHQLLQDESSYIFVSVTQEAEREEFFDETRRLCDLRLFQPFLKVIEPVGNREEKILNREIGFAIGMPVCEFDMVKDPEVQDFRRNILNVCKEAVDLRDLNSPHSRAMYVYPPNVESSPELPKHIYNKLDKGQIIVVIWVIVSPNNDKQKYTLKINHDCVPEQVIAEAIRKKTRSMLLSSEQLKLCVLEYQGKYILKVCGCDEYFLEKYPLSQYKYIRSCIMLGRMPNLMLMAKESLYSQLPMDCFTMPSYSRRISTATPYMNGETSTKSLWVINSALRIKILCATYVNVNIRDIDKIYVRTGIYHGGEPLCDNVNTQRVPCSNPRWNEWLNYDIYIPDLPRAARLCLSICSVKGRKGAKEEHCPLAWGNINLFDYTDTLVSGKMALNLWPVPHGLEDLLNPIGVTGSNPNKETPCLELEFDWFSSVVKFPDMSVIEEHANWSVSREAGFSYSHAGLSNRLARDNELRENDKEQLKAISTRDPLSEITEQEKDFLWSHRHYCVTIPEILPKLLLSVKWNSRDEVAQMYCLVKDWPPIKPEQAMELLDCNYPDPMVRGFAVRCLEKYLTDDKLSQYLIQLVQVLKYEQYLDNLLVRFLLKKALTNQRIGHFFFWHLKSEMHNKTVSQRFGLLLESYCRACGMYLKHLNRQVEAMEKLINLTDILKQEKKDETQKVQMKFLVEQMRRPDFMDALQGFLSPLNPAHQLGNLRLEECRIMSSAKRPLWLNWENPDIMSELLFQNNEIIFKNGDDLRQDMLTLQIIRIMENIWQNQGLDLRMLPYGCLSIGDCVGLIEVVRNSHTIMQIQCKGGLKGALQFNSHTLHQWLKDKNKGEIYDAAIDLFTRSCAGYCVATFILGIGDRHNSNIMVKDDGQLFHIDFGHFLDHKKKKFGYKRERVPFVLTQDFLIVISKGAQECTKTREFERFQEMCYKAYLAIRQHANLFINLFSMMLGSGMPELQSFDDIAYIRKTLALDKTEQEALEYFMKQMNDAHHGGWTTKMDWIFHTIKQHALN'


muts0 gives the wild fastas for each of the mutations.


At the time of writing NetSurfp1.1 can cope with 10,000 characters per file. Where more are needed a number of consecutive files are created
labelled muts0, muts1, muts2 etc.


The upper limit is expected to increase with the advent of NetSurfP1.2 and I will endeavour to update this code at that time.
Errors are saved to an errors.txt file in the output directory.


""")

parser.add_argument("path", help = """
Path to a text file of mutations, one on each line. 
Each mutation should be of the form PIK3CA_E545K

The first part is the name of the gene. This must be one of the following:
Gene_name, Entrez_gene_ID, Gene_stable_ID,Protein_stable_ID,
RefSeq_peptide_ID, Transcript_stable_ID,UniProtKB/Swiss-Prot_ID,
UniProtKB/TrEMBL_ID.

The next three terms refer to the name of the wild amino acid, the position of the amino-
acid and the name of the mutated amino acid.

""")
parser.add_argument("gene_type", help = """
One of the following:
Gene_name, Entrez_gene_ID, Gene_stable_ID,Protein_stable_ID,
RefSeq_peptide_ID, Transcript_stable_ID,UniProtKB/Swiss-Prot_ID,
UniProtKB/TrEMBL_ID.
Please do not mix gene name types
""")
parser.add_argument("output_dire", help = """Directory to save
the output in.""")
args = parser.parse_args()
path = args.path
gene_type = ' '.join(args.gene_type.split('_'))
output_dire = args.output_dire

gene_type = ' '.join(args.gene_type.split('_'))
output_dire = args.output_dire

data_dire = os.path.dirname(os.path.abspath(__file__))

def get_multi_converter():
    '''multi_converter enables conversions between 'Entrez gene ID', 'Fasta', 'Gene name', 
    'Gene stable ID','Protein stable ID', 'RefSeq peptide ID', 'Transcript stable ID',
    'UniProtKB/Swiss-Prot ID', 'UniProtKB/TrEMBL ID'. It is constructed from information 
    downloaded from biomart and from David https://david.ncifcrf.gov/conversion.jsp
    '''
    return pd.DataFrame.from_csv(os.path.join(data_dire,'/its/home/skw24/mutation_converter/mutation_converter/gene_multi_converter.csv'))

def converter(gene,labela,labelb):
    """Enables conversion of genes/proteins between different biomart labels. 
        
    Input 

    gene - string identifier
    labela = type of identifier - should be in list 'Entrez gene ID', 'Fasta', 'Gene name', 'Gene stable ID',
       'Protein stable ID', 'RefSeq peptide ID', 'Transcript stable ID',
       'UniProtKB/Swiss-Prot ID', 'UniProtKB/TrEMBL ID'
    labelb = type of identifier to translate to - again should be in the same list

    
    Ouput

    returns np.nan if there is no available translation

    """
    try:
        L=list(set(multi_converter[multi_converter[labela]==gene][labelb])-set([np.nan]))

        if len(L)>0:
            return L
        else:
            return np.nan
    except IndexError:
        print(gene)
        return np.nan

class Mutation:

    multi_converter = pd.DataFrame.from_csv(os.path.join(data_dire,'gene_multi_converter.csv'))
    
    def converter(gene,labela,labelb):
        """Enables conversion of genes/proteins between different biomart labels. 

        Input 

        gene - string identifier
        labela = type of identifier - should be in list 'Entrez gene ID', 'Fasta', 'Gene name', 'Gene stable ID',
           'Protein stable ID', 'RefSeq peptide ID', 'Transcript stable ID',
           'UniProtKB/Swiss-Prot ID', 'UniProtKB/TrEMBL ID'
        labelb = type of identifier to translate to - again should be in the same list


        Ouput

        returns np.nan if there is no available translation

        """
        try:
            L=list(set(Mutation.multi_converter[
                Mutation.multi_converter[labela]==gene][labelb])-set([np.nan]))

            if len(L)>0:
                return L
            else:
                return np.nan
        except IndexError:
            print(gene)
            return np.nan
    
    no_fasta =   """No fasta has been found. Is the gene_type correct?
            Choose between 'Entrez gene ID', 'Gene name', 'Gene stable ID',
       'Protein stable ID', 'RefSeq peptide ID', 'Transcript stable ID',
       'UniProtKB/Swiss-Prot ID', 'UniProtKB/TrEMBL ID. 
       The database is not complete.
       If you can't find your gene, try using the Transcript stable ID instead
       """
    
    no_appropriate_fasta = """No appropriate fasta has been found. 
            This may be because none of the identified fastas are long enough
            or because they have a different wild amino-acid at the position given
            """

    def __init__(self, m, gene_type):
        
        self.mut = m
        mut = self.mut.split('_')
        self.gene = mut[0]
        self.wild = mut[1][0]
        self.mutated = mut[1][-1]
        self.position = int(mut[1][1:-1])-1
        fastas = Mutation.converter(self.gene,gene_type,'Fasta')
        if type(fastas)==float:#ie fastas == np.nan
            self.wild_fasta = ''
            self.mut_fasta = ''
            self.message = Mutation.no_fasta
          
        else:
            fastas1 = [i for i in fastas if len(i)> self.position]
            fastas2 = [i for i in fastas1 if i[self.position]==self.wild]
            
            if len(fastas2)>0:
                lengths = [len(i) for i in fastas2]
                self.wild_fasta = fastas2[lengths.index(min(lengths))]
                self.mut_fasta = self.wild_fasta[:self.position]+\
                self.mutated+self.wild_fasta[self.position+1:]  
                self.message = ''
            else:
                self.wild_fasta = ''
                self.mut_fasta = ''
                self.message = Mutation.no_appropriate_fasta
                

            

def make_multi_list(name,mylist,output_dire):
    '''Takes a list of netsurfp inputs and splits them into appropriate lengths'''
    n=0
    a = [i for i in mylist]
    while a:

        short_list = []
        length = len(a[-1])

        while length<99999 and a:
            short_list.append(a.pop())
            if a:
                length+=len(a[-1])
        with open(os.path.join(output_dire,name+str(n)),'w') as f:
            f.write('\n'.join(short_list))
        n+=1


def main():

    with open(path,'r') as f:

        try:
            mutation_list = f.read().split('\n')

        except FileNotFoundError:
            print('file not found')
        
        #take out duplicates
        mutation_list =list(pd.Series(mutation_list).drop_duplicates())        
        muts = [Mutation(i,gene_type) for i in mutation_list]
        wilds = [i.wild_fasta for i in muts]
        mutants = [i.mut_fasta for i in muts]
        messages = [i.message for i in muts]

    #group all the wild fastas together so that the same one is not repeated twice
    wilds0 = pd.DataFrame(pd.Series(wilds,index = mutation_list),columns = ['Fastas'])
    wild_series = pd.Series(wilds0.groupby(by = 'Fastas').groups)
    wild_series1=pd.Series(wild_series.index)
    wild_series1.index = wild_series.map(lambda x:'|'.join(x))
    s2 = wild_series1.loc[wild_series1!='']
    wild_fasta_list = list('>'+s2.index+'\n'+s2)
    mutant_series = pd.Series(mutants,index = mutation_list)
    mutant_series = mutant_series[mutant_series!='']
    mut_fasta_list = list('>'+mutant_series.index+'\n'+mutant_series)

    #explain what has gone wrong with the process
    messages_series = pd.Series(messages, index = mutation_list)
    wrong = messages_series[messages_series!='']
    with open(os.path.join(output_dire,'errors.txt'),'w') as f:
        f.write('\n'.join(list(wrong.index+'\n'+wrong)))

    #split up the fasta lists into manageable lengths
    make_multi_list('wild',wild_fasta_list,output_dire)
    make_multi_list('mutant',mut_fasta_list,output_dire)

if '__name__' == '__main__':
    main()

