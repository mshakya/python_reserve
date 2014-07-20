#! /usr/bin/env python

"""
contains function to parse genbank .gbk and also .gbs file format
"""

from Bio import SeqIO
from Bio import SeqFeature
import re
from os import walk


class GenbankPro():
    """
    classes to parse genbank functions
    """
    def __init__(self, genbank_file):
        """create a new instance of a command line wrapper object."""
        self.read = SeqIO.parse(open(genbank_file, 'r'), 'genbank')

    def get16S(self, fasta_file):
        """ function to extract 16S rRNA genes from genbank file of a complete ncbi genome
        """
        fasta_handle = open(fasta_file, 'w')
        for genome in self.read:
            for gene in genome.features:
                if gene.type == 'source':
                    if gene.qualifiers['db_xref'][0].split(':')[0] == 'taxon':
                        taxonid = gene.qualifiers['db_xref'][0].split(':')[1]
                if gene.type == "rRNA":
                    if 'product' in gene.qualifiers:
                        if '16S' in gene.qualifiers['product'][0]:
                            start = gene.location.nofuzzy_start
                            end = gene.location.nofuzzy_end
                            if 'db_xref' in gene.qualifiers:
                                gi = []
                                gi = str(gene.qualifiers['db_xref'])
                                gi = gi.split(":")[1]
                                gi = gi.split("'")[0]
                                fasta_handle.write(">gi|%s|16S|taxid|%s|%s|%s\n%s\n" % (gi, taxonid, genome.description, genome.id, genome.seq[start:end]))
                            else:
                                fasta_handle.write(">gi|NoGenID|16S|taxid|%s|%s|%s\n%s\n" % (taxonid, genome.description, genome.id, genome.seq[start:end]))
        fasta_handle.close()

    def get16S2(self, fasta_file, out_fasta):
        """
        function to extract 16S rRNA genes of organisms housing CDS listed in given fasta_file
        """

        fasta_handle = open('out_fasta', 'w')
        gi_list = [str(seq.id.split('|')[1]) for seq in SeqIO.parse(open(fasta_file, 'r'), 'fasta')]

        genome_list = []
        for genome in self.read:
            for gene in genome.features:
                if gene.type == 'CDS':
                    if gene.qualifiers['db_xref'][0].split(':')[0] == 'GI':
                        if gene.qualifiers['db_xref'][0].split(':')[1] in gi_list:
                            genome_list.append(genome.id)

        for genome in self.read:
            if genome.id in genome_list:
                for gene in genome.features:
                    print gene
                    if gene.type == "rRNA":
                        if 'product' in gene.qualifiers:
                            if '16S' in gene.qualifiers['product'][0]:
                                start = gene.location.nofuzzy_start
                                end = gene.location.nofuzzy_end
                                if 'db_xref' in gene.qualifiers:
                                    gi = []
                                    gi = str(gene.qualifiers['db_xref'])
                                    gi = gi.split(":")[1]
                                    gi = gi.split("'")[0]
                                    fasta_handle.write(">gi|16S|taxid|%s|%s|%s\n%s\n" % (gi, genome.description, genome.id, genome.seq[start:end]))
                                else:
                                    fasta_handle.write(">gi|NoGenID|16S|taxid|%s|%s\n%s\n" % (taxonid, genome.description, genome.id, genome.seq[start:end]))

        fasta_handle.close()

    def get_cds(self, fasta_file):
        """
        function to extract all CDS from given genbank file
        """
        output_handle = open(fasta_file, 'w')
        for genome in self.read:
            for gene in genome.features:
                if gene.type == 'CDS':
                    assert len(gene.qualifiers['translation']) == 1
                    output_handle.write(">gi|%s|ref|%s| %s [%s]\n%s\n" % (
                            gene.qualifiers['db_xref'][0].split(':')[1],
                            gene.qualifiers['protein_id'][0],
                            gene.qualifiers['product'][0],
                            genome.description,
                            gene.qualifiers['translation'][0]))

        output_handle.close()


    def taxid2bioproj(self, taxid):
        """
        the function outputs a list of bioproject id for list of taxid provided from a given .gbs file
        """

        gbans = self.read
        bproj_list = []
        for genome in gbans:
            for feature in genome.features:
                if feature.qualifiers['db_xref'][0].split(":")[0] == 'taxon':
                    taxid_gbs = feature.qualifiers['db_xref'][0].split(":")[1]
                elif feature.qualifiers['db_xref'][1].split(":")[0] == 'taxon':
                    taxid_gbs = feature.qualifiers['db_xref'][0].split(":")[1]
                else:
                    print "there is no taxon specified in the first two db_xref of %s" % genome.source
                if taxid_gbs in taxid:
                    bproj = genome.dbxrefs[0].split(":")[1]
                    bproj_no = re.findall(r'\d+', bproj)
                    bproj_list.append(bproj_no[0])
        return bproj_list

    def taxid2folder(self, taxid):
        """
        creates a list of folder names that are similar to that listed in the ncbi ftp site

        Example:Acetobacter_pasteurianus_IFO_3283_26_uid158531/
        """
        gbans = self.read
        folder_list = []
        for genome in gbans:
            #organism name
            source_name = genome.annotations['organism']
            for feature in genome.features:
                if feature.qualifiers['db_xref'][0].split(":")[0] == 'taxon':
                    taxid_gbs = feature.qualifiers['db_xref'][0].split(":")[1]
                elif feature.qualifiers['db_xref'][1].split(":")[0] == 'taxon':
                    taxid_gbs = feature.qualifiers['db_xref'][1].split(":")[1]
                else:
                    print "there is no taxon specified in the first two db_xref of %s" % genome.source

            if taxid_gbs in taxid:
                #getting the uid number at the end of the folder name
                bproj = genome.dbxrefs[0].split(":")[1]
                bproj_no = re.findall(r'\d+', bproj)

                #to create folder name by substitutin - and dots and sp.

                #this if statement was created for some genome name that had = to sign at the end
                if source_name.find('=') != -1:
                    source_name = source_name.split(' =')[0]

                #for replacing sp.
                sname_replaced_v1 = source_name.replace(' sp.', '')
                #for replacing str.
                sname_replaced_v2 = sname_replaced_v1.replace(' str.', '')
                #for replacing St.
                sname_replaced_v3 = sname_replaced_v2.replace(' St.', '')

                # this is a special case for Beijerinckia_indica_subsp__indica_ATCC_9039_uid59057, need to find a way to fix this
                sname_replaced_v4 = sname_replaced_v3.replace(' subsp. indica', '')
                #special case due to folder name different thatn the genbank name
                #need to figure out how to use re.sub for this
                #sname_replaced_v3=re.sub(r'[' sp.',' str.',' St.']','',source_name)
                sname_replaced_v5 = re.sub(r'[\s,\-,\.,\#,\',\/]', '_', sname_replaced_v4)
                sname_replaced_v6 = sname_replaced_v5.replace('_DSM_122', '')
                sname_replaced = sname_replaced_v6.replace('extorquens_CM4', 'chloromethanicum_CM4')
                #sname_replaced=re.sub(r'[\s,\-,\.]','_',source_name)
                folder_name = sname_replaced+'_uid'+bproj_no[0]+'/'
                folder_list.append(folder_name)

        return folder_list

    def taxid2folder2(self, folder, taxid):

        """
        creates a list of folder names that are similar to that listed in the ncbi ftp site

        Example:Acetobacter_pasteurianus_IFO_3283_26_uid158531/
        """
    #read in the genbank file
        gbans = self.read
        folder_list = []
        #for creating a list of genome folders that were downloaded
        genome_list = []
        for (dirpath, dirnames, filenames) in walk(folder):
            genome_list.extend(dirnames)

        #for creating a list of bioproject based on the taxid
        for genome in gbans:
            for feature in genome.features:
                if feature.qualifiers['db_xref'][0].split(":")[0] == 'taxon':
                    taxid_gbs = feature.qualifiers['db_xref'][0].split(":")[1]
                elif feature.qualifiers['db_xref'][1].split(":")[0] == 'taxon':
                    taxid_gbs = feature.qualifiers['db_xref'][1].split(":")[1]
                else:
                    print "there is no taxon specified in the first two db_xref of %s" % genome.source

            if taxid_gbs in taxid:
                #getting the uid number at the end of the folder name
                bproj = genome.dbxrefs[0].split(":")[1]
                bproj_no = re.findall(r'\d+', bproj)
                bproj_uid = 'uid'+str(bproj_no[0])

                for genome_folders in genome_list:
                    if bproj_uid in genome_folders:
                        folder_list.append(genome_folders)

        return folder_list

    def gbankconv(self, gbank_file, gbank_type):
        """
        function to extract features from genbank file

        For example: gbank_file can be a concatenated genbank file,
        gbank_type=type to be extracted, for e.g. source or CDS

        """

        gbank = self.read

        out_file = gbank_file.split('.')[0]+'_%s.gbk' % gbank_type
        out_handle = open(out_file, 'w+')
        for genome in gbank:
            genome.features = [f for f in genome.features if f.type == gbank_type]
            SeqIO.write(genome, out_handle, "genbank")
        out_handle.close()

    def acceNgi(self, gbk_file, output_file):
        """
        function to read in genbank file and produce two column file, gi and accesssion
        """
        all_genomes = self.read

        out_handle = open(output_file, 'w')
        for genomes in all_genomes:
            genomes_acc = genomes.id
            for sequences in genomes.features:
                if sequences.type == 'CDS':
                    for records in sequences.qualifiers['db_xref']:
                        if records.split(':')[0] == 'GI':
                            gi_number = records.split(':')[1]
                            out_handle.write('%s\t%s\n' % (genomes_acc, gi_number))

        out_handle.close()

    def gbankTaxa(self):
        """
        extract taxonomic information from genbank file
        returns list of tuples that contains a list of taxonomy and name of organism
        """
        organism_list = []
        taxa_list = []
        gbank = self.read
        for genome in gbank:
            if genome.annotations['organism'] not in organism_list:
                organism_list.append(genome.annotations['organism'])
                taxa = genome.annotations['taxonomy'][3:4]
                taxa.append(genome.annotations['organism'])
                taxa_list.append(taxa)
        return taxa_list


def get16S(genbank_file, fasta_file):
    """ function to extract 16S rRNA genes from genbank file of a complete ncbi genome

    """
    from Bio import SeqIO, SeqFeature

    fasta_handle = open(fasta_file, 'w')
    gbank = SeqIO.parse(open(genbank_file, "rU"), "genbank")
    for genome in gbank:
        for gene in genome.features:
            if gene.type == 'source':
                if gene.qualifiers['db_xref'][0].split(':')[0] == 'taxon':
                    taxonid = gene.qualifiers['db_xref'][0].split(':')[1]
            if gene.type == "rRNA":
                if 'product' in gene.qualifiers:
                    if '16S' in gene.qualifiers['product'][0]:
                        start = gene.location.nofuzzy_start
                        end = gene.location.nofuzzy_end
                        if 'db_xref' in gene.qualifiers:
                            gi = []
                            gi = str(gene.qualifiers['db_xref'])
                            gi = gi.split(":")[1]
                            gi = gi.split("'")[0]
                            fasta_handle.write(">gi|%s|16S|taxid|%s|%s|%s\n%s\n" % (gi, taxonid, genome.description, genome.id, genome.seq[start:end]))
                        else:
                            fasta_handle.write(">gi|NoGenID|16S|taxid|%s|%s|%s\n%s\n" % (taxonid, genome.description, genome.id, genome.seq[start:end]))
    fasta_handle.close()


def taxid2bioproj(gbs_file, taxid):
    """
    the function outputs a list of bioproject id for list of taxid provided from a given .gbs file
    """
    import re
    from Bio import SeqIO

    gbans = SeqIO.parse(open(gbs_file, 'r'), 'genbank')
    bproj_list = []
    for genome in gbans:
        for feature in genome.features:
            if feature.qualifiers['db_xref'][0].split(":")[0] == 'taxon':
                taxid_gbs = feature.qualifiers['db_xref'][0].split(":")[1]
            elif feature.qualifiers['db_xref'][1].split(":")[0] == 'taxon':
                taxid_gbs = feature.qualifiers['db_xref'][0].split(":")[1]
            else:
                print "there is no taxon specified in the first two db_xref of %s" % genome.source
            if taxid_gbs in taxid:
                bproj = genome.dbxrefs[0].split(":")[1]
                bproj_no = re.findall(r'\d+', bproj)
                bproj_list.append(bproj_no[0])
    return bproj_list


def taxid2folder(gbs_file, taxid):
    """
    creates a list of folder names that are similar to that listed in the ncbi ftp site

    Example:Acetobacter_pasteurianus_IFO_3283_26_uid158531/
    """
    from Bio import SeqIO
    import re

    gbans = SeqIO.parse(open(gbs_file, 'r'), 'genbank')
    folder_list = []
    for genome in gbans:
        #organism name
        source_name = genome.annotations['organism']
        #print source_name
        for feature in genome.features:
            if feature.qualifiers['db_xref'][0].split(":")[0] == 'taxon':
                taxid_gbs = feature.qualifiers['db_xref'][0].split(":")[1]
            elif feature.qualifiers['db_xref'][1].split(":")[0] == 'taxon':
                taxid_gbs = feature.qualifiers['db_xref'][1].split(":")[1]
            else:
                print "there is no taxon specified in the first two db_xref of %s" % genome.source

        if taxid_gbs in taxid:
            #getting the uid number at the end of the folder name
            bproj = genome.dbxrefs[0].split(":")[1]
            bproj_no = re.findall(r'\d+', bproj)

            #to create folder name by substitutin - and dots and sp.

            #this if statement was created for some genome name that had = to sign at the end
            if source_name.find('=') != -1:
                source_name = source_name.split(' =')[0]

            #for replacing sp.
            sname_replaced_v1 = source_name.replace(' sp.', '')
            #for replacing str.
            sname_replaced_v2 = sname_replaced_v1.replace(' str.', '')
            #for replacing St.
            sname_replaced_v3 = sname_replaced_v2.replace(' St.', '')

            # this is a special case for Beijerinckia_indica_subsp__indica_ATCC_9039_uid59057, need to find a way to fix this
            sname_replaced_v4 = sname_replaced_v3.replace(' subsp. indica', '')
            #special case due to folder name different thatn the genbank name

            #need to figure out how to use re.sub for this
            #sname_replaced_v3=re.sub(r'[' sp.',' str.',' St.']','',source_name)
            sname_replaced_v5 = re.sub(r'[\s,\-,\.,\#,\',\/]', '_', sname_replaced_v4)
            sname_replaced_v6 = sname_replaced_v5.replace('_DSM_122', '')
            sname_replaced = sname_replaced_v6.replace('extorquens_CM4', 'chloromethanicum_CM4')
            #sname_replaced=re.sub(r'[\s,\-,\.]','_',source_name)
            folder_name = sname_replaced+'_uid'+bproj_no[0]+'/'
            folder_list.append(folder_name)

    return folder_list


def taxid2folder2(gbs_file, folder, taxid):

    """
    creates a list of folder names that are similar to that listed in the ncbi ftp site

    Example:Acetobacter_pasteurianus_IFO_3283_26_uid158531/
    """

    from Bio import SeqIO
    from os import walk
    import re

    #read in the genbank file
    gbans = SeqIO.parse(open(gbs_file, 'r'), 'genbank')
    folder_list = []

    #for creating a list of genome folders that were downloaded
    genome_list = []
    for (dirpath, dirnames, filenames) in walk(folder):
        genome_list.extend(dirnames)

    #for creating a list of bioproject based on the taxid
    for genome in gbans:
        for feature in genome.features:
            if feature.qualifiers['db_xref'][0].split(":")[0] == 'taxon':
                taxid_gbs = feature.qualifiers['db_xref'][0].split(":")[1]
            elif feature.qualifiers['db_xref'][1].split(":")[0] == 'taxon':
                taxid_gbs = feature.qualifiers['db_xref'][1].split(":")[1]
            else:
                print "there is no taxon specified in the first two db_xref of %s" % genome.source

        if taxid_gbs in taxid:
            #getting the uid number at the end of the folder name
            bproj = genome.dbxrefs[0].split(":")[1]
            bproj_no = re.findall(r'\d+', bproj)
            bproj_uid = 'uid'+str(bproj_no[0])

            for genome_folders in genome_list:
                if bproj_uid in genome_folders:
                    folder_list.append(genome_folders)

    return folder_list


def gbankconv(gbank_file, gbank_type):
    """
    function to subset genbank file

    For example: gbank_file can be a concatenated genbank file, gbank_type=type to be extracted, for e.g. source or CDS

    """

    from Bio import SeqIO
    gbank = SeqIO.parse(open(gbank_file, "rU"), "genbank")

    out_file = gbank_file.split('.')[0]+'_%s.gbk' % gbank_type
    out_handle = open(out_file, 'w+')
    for genome in gbank:
        genome.features = [f for f in genome.features if f.type == gbank_type]
        SeqIO.write(genome, out_handle, "genbank")
    out_handle.close()


def acceNgi(gbk_file, output_file):
    """
    function to read in genbank file and produce two column file, one with gi and accesssion
    """
    from Bio import SeqIO
    all_genomes = SeqIO.parse(gbk_file, 'genbank')

    out_handle = open(output_file, 'w')
    for genomes in all_genomes:
        genomes_acc = genomes.id
        for sequences in genomes.features:
            if sequences.type == 'CDS':
                for records in sequences.qualifiers['db_xref']:
                    if records.split(':')[0] == 'GI':
                        gi_number = records.split(':')[1]
                        out_handle.write('%s\t%s\n' % (genomes_acc, gi_number))

    out_handle.close()


def gbankTaxa(gbank_file):
    """
    extract taxonomic information from genbank file
    returns list of tuples that contains a list of taxonomy and name of organism
    """
    from Bio import SeqIO
    organism_list = []
    taxa_list = []
    gbank = SeqIO.parse(open(gbank_file, "rU"), "genbank")
    for genome in gbank:
        if genome.annotations['organism'] not in organism_list:
            organism_list.append(genome.annotations['organism'])
            taxa = genome.annotations['taxonomy'][3:4]
            taxa.append(genome.annotations['organism'])
            taxa_list.append(taxa)
    return taxa_list

