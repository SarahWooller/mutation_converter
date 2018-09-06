import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="mutation_converter",
    version="0.0.1",
    author="Sarah Wooller",
    author_email="s.k.wooller@sussex.ac.uk",
    description="Conversion of mutations into wild and mutated fastas",
    long_description="""The NetSurfP server predicts the surface accessibility and secondary structure of amino acids 
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
""",
    long_description_content_type="text/markdown",
    url="https://github.com/pypa/sampleproject",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
