#!/usr/bin/env python

#File created on Feb 16 2014
__author__ = "Migun Shakya"
__email__ = "microbeatic@gmail.com"
__version__ = "0.1"
__license__ = "The MIT License (MIT)"


"""
list of functions that deals with fasta
"""


#--- standard library imports
#


#third party packages
from Bio import SeqIO


class FastaPro():
    """
    class with methods to work with fasta files
    """
    def __init__(self, fasta_file):
        """create a new instance of a command line wrapper object."""
        self.read = SeqIO.parse(open(fasta_file, 'r'), 'fasta')

    def alpha_lysolytic(self, lysogen_fasta, lytic_fasta):
        '''
        function that distinguishes the alphaproteobacteria viruses as
        lysogentic and lytic phage
        '''
        alpha_virus = self.read
        lysogen_handle = open(lysogen_fasta, 'w')
        lytic_handle = open(lytic_fasta, 'w')

        lysogenic_list = ['Azospirillum phage Cd', 'Liberibacter phage SC1',
                          'Liberibacter phage SC2', 'Rhizobium phage 16-3',
                          'Rhizobium phage RR1-A', 'Rhizobium phage RR1-B',
                          'Rhodobacter phage RC1', 'Rhodobacter phage RcapMu',
                          'Rhodobacter phage RcapNL', 'Rhodovulum phage RS1 ']
        for seq in alpha_virus:
            if seq.description.split('[')[1].split(']')[0] in lysogenic_list:
                lysogen_handle.write('>'+seq.description+'\n'+seq.seq)
            else:
                lytic_handle.write('>'+seq.description+'\n'+seq.seq)
        lytic_handle.close()
        lysogen_handle.close()

    def extract_fasta_panseq(self, list_file, out_fasta):
        '''
        function to extract sequences listed in a text file from a fasta_file
        '''

        source_seqs = self.read
        id_list = []
        with open(list_file, 'r') as file_list:
            for line in file_list:
                id_list.append(line.split('\n')[0])
        with open(out_fasta, 'w') as out_file:
            for seq in source_seqs:
                if seq.id.split('|')[1] in id_list:
                    out_file.write('>'+str(seq.description)+'\n')
                    out_file.write(str(seq.seq)+'\n')


def seqlist(fasta_file, attribute_file):
    """
    extracts sequence description that matches strings in the attribute_file and gives out a list

    fasta_file: file from which the sequence header that contains strings in attribute_list

    attribute_file: file containing column of strings to be searched
    """
    seq_descs = []

    for attributes in open(attribute_file, 'r'):
        in_fasta = SeqIO.parse(fasta_file, 'fasta')
        for sequences in in_fasta:
            if attributes in sequences.description:
                seq_descs.append(sequences.description)
    return seq_descs


def removeseq(fasta_file, seq_list):
    """
    create '*_v1.fasta' after removing sequence from fasta_file with header that matches the one provided in seq_list.

    fasta_file=fasta file from which the sequence whose header matches that in the seq_list are removed
    seq_list: list of sequence header to be removed from the fasta_file
    """

    out_fasta = open(fasta_file+'_v1.fasta','w')
    in_fasta=SeqIO.parse(fasta_file,'fasta')
    for sequences in in_fasta:
        if sequences.description not in seq_list:
            out_fasta.write('>'+str(sequences.description)+'\n')
            out_fasta.write(str(sequences.seq)+'\n')
    out_fasta.close()


def fastadict(fasta_file, separ, pos):
    """
    creates dictionary with a part of sequence header as key and sequence as value

    fasta_file=files to be read
    sep=separation that is used in the sequence header
    pos=position of header to be used as key based on sep separation

    For example fastadict(fasta_file,'|',3)  will create a dictionary with key as
    the segement from 4th position
    """
    fasta_dict = {}
    fasta_handle = SeqIO.parse(fasta_file, 'fasta')
    for sequence in fasta_handle:
        fasta_dict[sequence.description.split(separ)[pos]] = sequence.seq
    return fasta_dict


def renameseq(fasta_file):
    """
    renames fasta sequences with series of number indices

    'sequence' was added to the numbers to work with boxshade function,
    """
    out_file = fasta_file.split('.')[0]
    out_fasta = open(out_file+'_renme.faa', 'w')
    fasta_handle = SeqIO.parse(fasta_file, 'fasta')
    i = 0
    for sequences in fasta_handle:
        out_fasta.write('>'+str(i)+'sequenc'+'\n')
        out_fasta.write(str(sequences.seq)+'\n')
        i = i+1
    out_fasta.close()
    return str(out_file+'_renme.faa')


def gicsv2fasta(csv_file, source_fasta, result_fasta):
    """ extracts fasta file given the csv file that contains gi number in the second column
    """
    gi_file = open(csv_file, 'r')
    out_fasta = open(result_fasta, 'w')
    gi_list = []
    for line in gi_file:
        gi_list.append(line.split(',')[1].split('\n')[0])
    open_fasta = SeqIO.parse(open(source_fasta), 'fasta')
    for seq in open_fasta:
        if seq.id.split('|')[1] in gi_list:
            SeqIO.write(seq, out_fasta, 'fasta')
    out_fasta.close()
    gi_file.close()


def duplicate(fasta_file):
    """inputs a fasta file, outputs fasta file with
    duplicate sequences removed but sequence description combined in the header
    """
    sequences = {}
    for seq_record in SeqIO.parse(fasta_file, 'fasta'):
                sequence = str(seq_record.seq).upper()
                if sequence not in sequences:
                        sequences[sequence] = seq_record.id
                else:
                        sequences[sequence] + = '_'+seq_record.id

    out_file = fasta_file.split('.')[0]+'_uniq.faa'
    output_file = open(out_file, "w+")
    #Just Read the Hash Table and write on the file as a fasta format
    for sequence in sequences:
                output_file.write(">"+sequences[sequence]+"\n"+sequence+"\n")
    output_file.close()


def list2fasta(seqlist, source_fastafile, out_fastafile):
    """extracts fasta based on the list provided
    """

    result_handle = open(out_fastafile, 'w')
    seqiter = SeqIO.parse(open(source_fastafile, 'r'), 'fasta')
    for seq in seqiter:
        if seq.id in seqlist:
            SeqIO.write(seq, result_handle, 'fasta')

    result_handle.close()


def extract_fasta(in_fasta, out_fasta, start, end):
    """function that extracts the region of sequence
    based on start an end position provided
    """
    from Bio import SeqIO

    fasta_handle = SeqIO.parse(in_fasta, 'fasta')
    out_handle = open(out_fasta, 'w')
    for sequence in fasta_handle:
        new_sequence = sequence[start:end]
        out_handle.write('>' + str(sequence.description) + '\n' + str(new_sequence) + '\n')
    out_handle.close()


def extract_fasta(in_fasta, acc_id, start, end):
    """function that extracts the region of sequence
    based on start an end position provided
    """
    from Bio import SeqIO

    fasta_handle = SeqIO.parse(in_fasta, 'fasta')
    #out_handle = open(out_fasta, 'w')
    sorted_start_end = sorted([start, end])
    for sequence in fasta_handle:
        if str(sequence.id.split('|')[3]) == str(acc_id):
            new_sequence = sequence.seq[sorted_start_end[0]:sorted_start_end[1]]
    return str(new_sequence), sequence.description



