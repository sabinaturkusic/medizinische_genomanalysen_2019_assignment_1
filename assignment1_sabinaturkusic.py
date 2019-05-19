import mysql.connector
from os import path
import pysam

__author__ = 'Sabina Turkusic'

##
## Concept:
## TODO
##


class Assignment1:

    #Constructor
    def __init__(self, gene, genome_reference, file_name, bam_file):
        self.gene = gene
        self.genome_reference = genome_reference
        self.file_name = file_name #ucsc_name
        self.bam_file = bam_file
        self.alignment_file = pysam.AlignmentFile(self.bam_file, "rb")


    def download_gene_coordinates(self):
        print("Connecting to UCSC to fetch data")
        
        ## Open connection
        cnx = mysql.connector.connect(host='genome-mysql.cse.ucsc.edu', user='genomep', passwd='password', db=self.genome_reference)
        
        ## Get cursor
        cursor = cnx.cursor()
        
        ## Build query fields
        query_fields = ["refGene.name2",
                        "refGene.name",
                        "refGene.chrom",
                        "refGene.txStart",
                        "refGene.txEnd",
                        "refGene.strand",
                        "refGene.exonCount",
                        "refGene.exonStarts",
                        "refGene.exonEnds"]
        
        ## Build query
        query = "SELECT DISTINCT %s from refGene" % ",".join(query_fields)
        
        ## Execute query
        cursor.execute(query)
        
        ## Write to file
        ## TODO this may need some work
        with open(self.file_name, "w") as fh:
            for row in cursor:
                fh.write(str(row) + "\n")
    
            
        ## Close cursor & connection
        cursor.close()
        cnx.close()
        
        print("Done fetching data")

        
    def get_coordinates_of_gene(self):
        ## Use UCSC file
        count = 0
        with open(self.file_name, "r") as UCSC_file:
            for line in UCSC_file:
                if self.gene in line:
                    count+=1
                    print(line.rstrip("\n"))
        print("There are", count, "entries for", self.gene, "and only first entry will be used.")

        #use only first entry to get gene coordinates
        with open(self.file_name, "r") as UCSC_file:
            for line in UCSC_file:
                if self.gene in line:
                    line_split = line.replace(")", "").replace("(", "").replace("'", "")
                    line_split = line_split.split(", ")
                    break

        gene_symbol = line_split[0:2]
        chromosome = line_split[2]
        start = int(line_split[3])
        end = int(line_split[4])
        exons = line_split[6]

        print("GENE COORDINATES OF FIRST ENTRY: Chromosome:", chromosome, ", Start-position:", start, ", Stop-position:", end)

        #Gene Symbol

        print("The gene symbol is: ", gene_symbol)

        #Exons

        print("The number of exons is: ", exons)

        #Reads

        reads = self.alignment_file.fetch(chromosome, start, end)

        #Prpperly Paired

        count_pair = 0
        count_reads = 0
        for read in reads:
            count_reads += 1
            if read.is_proper_pair:
                count_pair += 1
        print("In total there are", count_reads, "reads and", count_pair, "are properly paired")


        #Reads with Indels
        reads = self.alignment_file.fetch(chromosome, start, end)
        indel_list = []
        for read in reads:
            if not read.is_unmapped:
                cig = read.cigartuples
                for (operation, length) in cig:
                    if (operation == 1) or (operation == 2):
                        indel_list.append(read)
        print("There are ", len(indel_list), "reads with indels.")


        #Number of mapped Reads
        reads = self.alignment_file.fetch(chromosome, start, end)
        count_mapped = 0
        for read in reads:
            if not read.is_unmapped:
                count_mapped += 1
        print("Number of mapped reads:", count_mapped)


        #Total Coverage - läuft länger
        #total_cov = []
        #for pileup in self.alignment_file.pileup(chromosome):
        #    total_cov.append(pileup.n)

        #average_total_coverage = sum(total_cov) / float(len(total_cov))
        #print("The average total coverage is: ", average_total_coverage)

        #Average Gene Coverage
        gene_cov = []
        for pileup in self.alignment_file.pileup(chromosome, start, end):
            gene_cov.append(pileup.n)

        average_gene_coverage = sum(gene_cov) / float(len(gene_cov))
        print("The average gene coverage is: ", average_gene_coverage)


    def get_sam_header(self):
        header_dict = self.alignment_file.header
        #print(header_dict)
        #print(samfile.header.keys())


def main():
    print("Assignment 1")
    assignment1 = Assignment1("NRIP1", "hg38", "UCSC_file.txt", "chr21.bam" )

    #check if file exists
    if (path.exists("UCSC_file.txt")):
        print("Processing gene coordinates...")
    else:
        assignment1.download_gene_coordinates()  #download gene coordinates

    assignment1.get_coordinates_of_gene()
    assignment1.get_sam_header()

    print("Done with assignment 1")
    
        
if __name__ == '__main__':
    main()
