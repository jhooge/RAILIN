'''
Created on Apr 25, 2012

@author: jhooge
'''
from Definitions import one2Three
from Residue import AminoAcid, PastaResidue
from Table import Table
from copy import deepcopy
import Definitions as DEF

class FastaRecord(object):
    
    def __init__(self, title, sequence, similarity=None):
        self.title = title ## title string
        self.sequence = sequence ## list of one letter codes
        self.similarity = similarity ## similarity float value
        
    def toAminoAcids(self, fn):
        """
        Generates Residue.AminoAcid objects from FASTA sequence
        
        @param fn: bmrb reference filename 
        @type fn: str
         
        @return: list of Residue.AminoAcids
        @rtype: list
        """
        fh = FileHandler()
        amino_acids = fh.read_preset(fn)
        seq = []
        for one_let in list(self.sequence):
            three_let = one2Three(one_let)
            for aa in amino_acids:
                if aa.three_let == three_let:
                    seq.append(deepcopy(aa))
        return seq

    def __str__(self):
        """
        String representation of FASTA Record
        """
        string = self.title + "\n" + self.sequence
        return string
        
class FileHandler(object):
    '''
    classdocs
    '''
    
    def __init__(self):
        '''
        Constructor
        '''
    
    def read_fasta_entry(self, fn):
        """
        The first line in a FASTA record has to be a the title 
        starting with ">"
            
        Examples:
        >Template
        >Sequence_0|0.68
        
        @raise TypeError: If fn is not a a file in FASTA format
        
        @param fn: filename of sequence file in FASTA format
        @type fn: string
        
        @return: a FASTARecord object filled with data from FASTA file
        @rtype: FASTARecord
        """
        infile = open(fn)
        line = infile.readline()
        
        if not line:
            # End of file
            return None
    
        # Double-check that it's a FASTA file by looking for the '>'.
        if not line.startswith(">"):
            raise TypeError("Not a FASTA file: %r" % line)
    
        # The title starts after the '>' and doesn't include any
        # trailing whitespace.
        title = line[1:].rstrip()
    
        # Read the sequence lines up to the blank line.
        sequence_lines = []
        while 1:
            # I know I'm at the end of the sequence lines when I reach the
            # blank line or the end of the file.  I can simplify the test
            # by calling rstrip and checking for the empty string.  If it's
            # a sequence line I'll still call rstrip to remove the final
            # newline and any trailing whitespace so I'll go ahead and
            # always rstring the line.
            line = infile.readline().rstrip()
            # remove whitespaces
            line = line.replace(' ', '') 
            
            if line == "":
                # Reached the end of the record or end of the file
                break
            sequence_lines.append(line)
        # Merge the lines together to make a single string containing the sequence
        # (This is faster than doing "sequence = sequence + line" if there are
        # more than a few lines)
        sequence = "".join(sequence_lines)
        record = FastaRecord(title, sequence)
        return record
    
    def read_pasta(self, pastafn, presetfn):
        import numpy as np
        '''
        Parses a pasta list and returns a list of accordingly
        filled Residue object. A Residue object is instanciated every time 
        a newline has been read.
        
        @attention: AT EOF THERE HAS TO BE A NEWLINE !
        
        @param pastafn: filename of residue list
        @type pastafn: string
        @param presetfn: filename of bmrb reference list
        @type presetfn: string
        
        @return: A list of Residue.PastaResidue filled with values of
        the residue list.
        @rtype: Residue.PastaResidue
        '''
        
        ## store preset atoms to check if res labels occur in preset
        preset = self.read_preset(presetfn)
        ref_atoms = []
        for a in preset:
            for atom in a.shifts.keys():
                if atom not in ref_atoms:
                    ref_atoms.append(atom)
                    
        infile = open(pastafn, 'r')
        content = infile.readlines()
        output = []
        
        header = content[0]
        lastLine = content[-1]
        if header != "#PASTA residue list\n" or content[1] != '\n':
            raise TypeError("\"%s\" is not a PASTA File!\n\
            First line should be \"#PASTA residue list\" followed by an empty line.\n\
            Given: %s" % (pastafn, header))
        if lastLine != '\n':
            raise TypeError("Last line has to be a carriage return!")
        
        for line in content[2:]:
            
                if(line != '\n' and not line.startswith('//')):
                    if(line.startswith('#')):
                        i = 0
                        info = line.strip().split('\t')
                    else:
                        tmp = line.strip().split()
                        info = info + tmp
                elif not line.startswith('//'):
                    shifts = info[4:]
                    pasta_no = int(info[0].split(':')[1].strip())
                    energy = int(float(info[1].split(':')[1].strip()))
                    seq_pos = int(info[2].split(':')[1].strip())
                    sec_struc = info[3].split(':')[1].strip()
                    shift_list = []

                    name = 'NAA'
                    name_im1 = 'NAA'

                    while(i < len(shifts)):
                        if shifts[i:i + 4][1].endswith('i-1'):
                            name_im1 = shifts[i:i + 4][0]
                        else:
                            name = shifts[i:i + 4][0]
                        shift_list.append(shifts[i:i + 4])
                        i = i + 4
                    
#                    aa = shift_list[0][0]
#                    print shift_list[0][3]
                    if (len(shift_list) > 0):
                        orig_no = shift_list[0][3]
                        orig_no = int(orig_no[1:-1])
                        
#                    
                    types = []
                    values = []
                    for entry in shift_list:
                        types.append(entry[1])
                        values.append(float(entry[2]))

                    shift_dict = dict(zip(types, values))
                    
                    types_i = []
                    types_im1 = []
                    values_i = []
                    values_im1 = []
                    
                    for entry in shift_list:
                        if entry[1].endswith('i-1'):
                            types_im1.append(entry[1])
                            values_im1.append(float(entry[2]))
                        else:
                            types_i.append(entry[1])
                            values_i.append(float(entry[2]))
                            
                    shifts_i = dict(zip(types_i, values_i))
                    shifts_im1 = dict(zip(types_im1, values_im1))
                    
                    ## check for invalid atom labels in pasta list
                    for atom in shifts_i.keys():
                        valid_atoms = ["C1", "C2", "C3", "C4", "C5", "C6","C7", "N15"]
                        if atom not in valid_atoms:
                            if atom not in ref_atoms:
                                raise TypeError("Residue no %i:\n \
                                Label \"%s\" not found in reference list.\n \
                                Please use IUPAC convention!" % (pasta_no, atom))
                    
                    for atom in shifts_im1.keys():
                        valid_atoms = ["C1", "C2", "C3", "C4", "C5", "C6", "C7", "N15"]
                        if atom[:-3] not in valid_atoms:
                            if atom[:-3] not in ref_atoms:
                                raise TypeError("Residue no %i:\n \
                                Label \"%s\" not found in reference list.\n \
                                Please use IUPAC convention!" % (pasta_no, atom))
                                
                    pasta_residue = PastaResidue()
                    if len(shift_list) > 0: 
                        pasta_residue.set_attribute('pasta_no', pasta_no)
                        pasta_residue.set_attribute('original_no', orig_no)
                        pasta_residue.set_attribute('energy', pasta_no)
                        pasta_residue.set_attribute('shifts_i', shifts_i)
                        pasta_residue.set_attribute('shifts_im1', shifts_im1)
                        pasta_residue.set_attribute('name', name)
                        pasta_residue.set_attribute('name_im1', name_im1)
                        pasta_residue.set_attribute('seq_pos', seq_pos)
                    output.append(pasta_residue)
                    
        for res in output:
            for key in res.shifts_i.keys():
                if res.shifts_i[key] < 0:
                    print 'Warning: In Residue %i, %s: %.2f will be ignored!' % (res.pasta_no, key, res.shifts_i[key])
                    del res.shifts_i[key]
            for key in res.shifts_im1.keys():
                if res.shifts_im1[key] < 0:
                    print 'Warning: In Residue %i, %s: %.2f will be ignored!' % (res.pasta_no, key, res.shifts_im1[key])
                    del res.shifts_im1[key]
        
#            if np.array([np.array(res.shifts_i.values())<0]).any()\
#            or np.array([np.array(res.shifts_im1.values())<0]).any():
#                raise FileFormatError('Negative chemical shift values found in!')
            if len(res.shifts_i.values()) + len(res.shifts_im1.values()) == 0:
                print 'Warning: Residue %i contains no measurements!' % res.pasta_no
#                print 'Warning: Residue %i will be removed!' % res.pasta_no
#                output.remove(res)
#                raise FileFormatError('Residue without chemical shift measurements found!')
        
        infile.close()
        return output
    
    def read_seq(self, seqfn, presetfn):
        """
        Reads a sequence file in FASTA format, populates fields 
        in Residue.AminoAcid object and returns a list of amino acids
        filled with reference values
        
        @param seqfn: filename for sequence file
        @type seqfn: string
        @param presetfn: filename for BMRB reference file
        @type presetfn: string
        @return: list of amino acid objects populated with reference
        shifts
        @rtype: list 
        """
        record = self.read_fasta_entry(seqfn)
        amino_acids = record.toAminoAcids(presetfn)
        return amino_acids
    
    def read_preset(self, fn):
        """
        Reads a preset file (e.g. BMRB) and returns it as  list of amino acid 
        objects with corresponding reference shifts
        
        @param fn: filename for BMRB reference file
        @type fn: string
        @return: list of amino acids with reference shifts
        @rtype: list
        """
        stats = self.to_stats(fn)
        amino_acids = self.statsToAminoAcids(stats)
        return amino_acids
    
    def read_Nhsqc(self, fn):
        """
        Reads Nhsqc file which is in table form and returns the
        first column which should be the identifiers (amino_acid type)
        for the data and the data in an numpy array.
        
        First line defines the header and all following lines
        the data values.
        
        @attention: First element in header needs to be called "Assignment"
        or the file is not recognized as Nhsqc file 
        
        Example:
        Assignment w1 w2 w3 w4 DataHeight
        ?-?-?-? 59.479 59.559 108.502 8.768 269785
        M1N-HN  59.543 30.310 120.345 7.819 225944
        ...
        
        @param fn: filename for Nhsqc file
        @type fn: string
        
        @return: tuple of identifier column and data array
        @rtype: tuple 
        """
        from numpy import array, zeros
        
        t = Table()
        t.read(fn)
        
        header = t[0]
        cols = len((t[2:][0]))
        rows = len(t[2:])
        
        if header[0] != "Assignment":
            raise TypeError("Not a NHSQC File")
        if len(header) != cols:
            raise TypeError("Header has more Identifiers than columns (whitespaces?)")
        
        identifiers = []
        data = zeros((rows, cols - 1))
        
        for i, row in enumerate(t[2:]):
            identifiers.append(row[0])
            data[i] = row[1:]
        
        return identifiers, data
    
    def statsToAminoAcids(self, stats):
        """
        Generates a list of amino acid objects from BMRB reference parameters.
        
        @param stats: dictionary of amino acid reference parameters.
                      keys: amino acid three letter code
                      values: list of mean and std values [mu, sigma]
        @type stats: dict
        
        @return: list of amino acid objects with respective parameters.
        @rtype: list
        """
        
        output = []
        for name in stats.keys():
            aa = AminoAcid()
            aa.set_attribute('one_let', DEF.three2One(name))
            aa.set_attribute('three_let', name)
            aa.set_attribute('shifts', stats[name])
            output.append(aa)
        
        # sort output according to Definitions file
        indexes = DEF.three_let2Index([aa.three_let for aa in output])
        tmp = output[:]
        
        i = 0
        for aa in output:
            tmp[indexes[i]] = aa
            i += 1
        
        output = tmp
        
        return output
    
    def write_pasta(self, fn, residues):
        """
        Generates a pasta fn with a given list of 
        Residue.PastaResidue objects
        
        @param fn: filename
        @type fn: string
        @param residues: list of Residue.PastaResidue objects
        @type residues: list
        """
        infile = open(fn, 'w')

        infile.write('#PASTA residue list\n\n')
        
        for res in residues:
            infile.write(str(res) + '\n')
        
        infile.close()
    
    def write_seq(self, fn, title, seq):
        """
        Writes sequence to file in FASTA file.
        
        @param fn: filename
        @type fn: string
        @param title: FASTA entry title
        @type title: string
        @param seq: list of Residue.AminoAcid objects
        @type seq: list
        """
        
        f = open(fn, 'w')
        f.write('>' + title + '\n')
        f.write(''.join([aa.one_let for aa in seq]))
        f.close()
    
    def write_fasta(self, fn, records, sequence):
        """
        Write a set of records plus template sequence to a fasta file
        
        @param fn: filename
        @type fn: string
        @param records: a list of FASTA Records
        @type records: list
        @param sequence: list of Residue.AminoAcid objects
        @type sequence: list
        """
        f = open(fn, 'w')
        string = self.__template_to_fasta(sequence)
        string += "\n"
        for record in records:
            string += str(record)
            string += "\n\n"
        f.write(string)
        f.close()
    
    def to_fasta_records(self, assignments, sequence, sopt_res):
        """
        Writes a set of residue to sequence assignments to a list 
        of FASTA records
        
        @param assignments: list of lists filled with tuples (Residue no, Sequence Position)
        @type assignments: list
        @param sequence: list of AminoAcid objects in sequence
        @type sequence: list
        @param sopt_res: list of PastaResidue objects in residue list
        @type sopt_res: list
        
        @return: records filled with assigned pasta records
        @rtype: list
        """
        records = []
        for i, assignment in enumerate(assignments):
            record = self.to_fasta_record(assignment, sequence, sopt_res, i)
            records.append(record)
            
        records.sort(key=lambda x: x.similarity, reverse=True)
        return records
    
    def to_fasta_record(self, assignment, sequence, residues, sequence_no=0):
        """
        Writes a given residue to sequence assignment to a FASTA record
        
        @param assignment: tuple of ([(Residue no, Sequence Position),...], similarity)
        @type assignment: tuple
        @param sequence: list of AminoAcid objects in sequence
        @type sequence: list
        @param residues: list of PastaResidue objects in residue list
        @type residues: list
        @param sequence_no: id for sequence title
        @type sequence_no: int
        
        @return: FASTA Record with title and sequence string derived 
        from assignment
        @rtype: FASTARecord
        """
        
        from Definitions import three2One
        line_break = 80
        assi = assignment[0]
        simi = assignment[1]
        title = '>Sequence_%i|%.2f' % (sequence_no, simi)
        seq = ['-'] * len(sequence)
        for i, j in assi:
            if i < len(residues):
                name = residues[i].name
                if len(name) > 3:
                    name = name[:3]
                if name != "NAA":
                    seq[j] = three2One(name)
        for k in range(line_break, len(sequence), line_break):
            seq.insert(k, '\n')
        seq = ''.join(seq)
        record = FastaRecord(title, seq, simi)
        return record
    
    def to_stats(self, fn):
        """
        Generates dictionary with keys and values from reference file
        (eg.: BMRB)
        Reference file has to have correct header:
        ['Res', 'Name', 'Atom', 'Count', 
        'Min.', 'Max.', 'Avg.', 'Std', 
        'Dev'] or it will not be recognized as reference file
        
        Keys in this dictionary are amino acid three letter codes
        and values are represented by another dictionary
        where keys are Atom types and there corresponding
        mean and standard deviation values in a list
        
        If you want to access the mean/standard deviation 
        value of the C-alpha atom of Cystein a call would look like this.
        
        mean = stats["CYS"]["CA"][0] 
        std = stats["CYS"]["CA"][1]
        
        @param fn: filename for BMRB reference file
        @type fn: string
        
        @return: Dictionary with amino acid shift values
        @rtype: dictionary
        """
        t = Table()
        t.read(fn)
        
        if t[0] != ['Res', 'Name', 'Atom', 'Count',
                    'Min.', 'Max.', 'Avg.', 'Std',
                    'Dev']:
            raise TypeError("Not a BMRB reference file: %r" % t[0])
        
        stats = {}
        
        for row in t:
        
            if not len(row) == 8: continue
        
            aa = row[0]
            if not aa in stats: stats[aa] = {}
            stats[aa][row[1]] = row[-2:]
            
        return stats
    
    def conf_to_pasta(self, fn, conf, sopt_res, amino_acids):
        """
        Generates a pasta file from given residue to sequence assignment
        
        @param fn: filename
        @type fn: string
        @param conf: Assignment tuples [(residue_no, seq_no)...]
        @type conf: list
        @param sopt_res: Residue Objects
        @type sopt_res: list
        @param amino_acids: Amino_acid objects 
        @attention: make sure to provide aa list with Cy* because it will be dealt with)
        @type amino_acids: list
        """
        from copy import deepcopy
        
        amino_acids = map(lambda x: x.one_let, amino_acids)
        k = amino_acids.index("B")
        del amino_acids[k]
        
        sorted_residues = deepcopy(sopt_res)
        assignments, similarity = conf
        for i, j in assignments:
            res = sorted_residues[i]
            res.set_attribute("seq_pos", j+1)
        
        sorted_residues.sort(key=lambda x: x.seq_pos)
        
        self.write_pasta(fn, sorted_residues)
    
    def __template_to_fasta(self, sequence):
        """
        Converts template sequence to FASTA formatted string
        
        @param sequence: list of AminoAcid objects in sequence
        @type sequence: list
        
        @return: FASTA formatted string
        @rtype: string
        """
        line_break = 80
        string = '>Template\n'
        seq = map(lambda x:x.one_let, sequence)
        for i in range(line_break, len(sequence), line_break):
            seq.insert(i, '\n')
        seq = ''.join(seq)
        string += seq + '\n'
        return string
    
    def assignment_matrix_from_assigned(self, pastafile, seqfile, statsfile):
        from numpy import zeros
        fh = FileHandler()
        sopt_res = fh.read_pasta(pastafile, statsfile)
        sequence = fh.read_seq(seqfile, statsfile)
        A = zeros((len(sopt_res), len(sequence)))
        for res in sopt_res:
            res_index = res.original_no - 1
            seq_index = res.seq_pos - 1
            if res.seq_pos > -1:
                A[res_index][seq_index] = 1
            else:
                print "No Assignment for ", res.original_no
        return A

    def expected_linking(self, residues):
        from numpy import zeros
        L = zeros((len(residues), len(residues)))
        for i, pred in enumerate(residues):
            for j, succ in enumerate(residues):
                if pred.seq_pos == succ.seq_pos - 1:
                    L[i][j] = 1
        return L

    def test_linking(self, A, L):
        for i, row in enumerate(L):
            if row.sum() == 0:
                print i, A[i].sum() == 0
    
    def read_bmrb_file(self, fn, bmrbfn):
        
        t = Table()
        t.read(bmrbfn)
        
        residues = []
        data = t.columns([5, 6, 7, 10])
        i=0
        j=0
        while j < len(data):
            pasta_no = data[i][0]
            name = data[i][1]
            shifts = {}
            while j < len(data) and data[j][0] == data[i][0]:
                nucleus = data[j][2]
                shift = data[j][3]
                if nucleus == "C":
                    nucleus = "CO"
                if nucleus == "N":
                    nucleus = "N15"
                if not (nucleus.startswith("H") or nucleus.startswith("NE")):
                    shifts[nucleus] = shift
                j+=1
            i = j
            res = PastaResidue()
            res.set_attribute("pasta_no", pasta_no)
            res.set_attribute("original_no", pasta_no)
            res.set_attribute("seq_pos", pasta_no)
            res.set_attribute("name", name)
            res.set_attribute("shifts_i", shifts)
            residues.append(res)
            
        ## set [i-1] shifts
        i=0
        j=1
        while j < len(residues):
            pred = residues[i]
            succ = residues[j]
            shifts_im1 = {}
            for key in pred.shifts_i.keys():
                if key != "N15":
                    shifts_im1[key+"i-1"] = pred.shifts_i[key]
            
            succ.set_attribute("shifts_im1",shifts_im1)
            succ.set_attribute("name_im1",pred.name)
            i+=1
            j+=1
        self.write_pasta(fn, residues)
    
    def add_linking_noise(self,residues):
        
        return residues
    
if __name__ == '__main__':
    from Pasta import Pasta
    from FileHandler import FileHandler
    from Definitions import three2One
    from numpy import array, dot, where
    from DataGenerator import DataGenerator
    from copy import deepcopy
    import pylab as pl
    import pickle
    
    names = ["SuperoxideDismutase"]
    for name in names:
        print name
        bmrbfn = "Datasets/%s/nmrstar/%s.str"%(name,name)
        pastafn = "Datasets/%s/residue_lists/%s.pasta"%(name,name)
        presetfn = "tests/reference_lists/bmrb.shift"
        seqfn = "Datasets/%s/sequences/%s.fasta"%(name,name)
        fh = FileHandler()
        dg = DataGenerator(presetfn)
        fh.read_bmrb_file("Datasets/%s/residue_lists/%s.pasta"%(name,name), bmrbfn)
        seq = fh.read_seq(seqfn, presetfn)
        residues = fh.read_pasta(pastafn, presetfn)
        res_reversed = sorted(residues, key=lambda x: x.pasta_no, reverse=True)
        for r in res_reversed:
            r.name = "NAA"
            r.name_im1 = "NAA"
            r.seq_pos = -1
        fh.write_pasta("Datasets/%s/residue_lists/%s_unassigned.pasta"%(name,name), res_reversed)
    
    
    
    
    
    
    
    
#    
#    optfile = "Datasets/Ubiquitin/Ub_opt_relabeled.pasta"
#    statsfile = 'tests/reference_lists/bmrb.shift'
#    seqfile = "Datasets/Ubiquitin/Ub.fasta"
#    
#    pasta = Pasta(optfile, statsfile, seqfile)
#    residues = pasta.residues
#    fh =  FileHandler()
#    
#    A = fh.assignment_matrix_from_assigned(optfile, seqfile, statsfile).astype("i")
##    L_expected = fh.linking_matrix_from_assigned(A).astype("i")
#    L_expected = fh.expected_linking(residues).astype("i")
#    L_estimated = pasta.linking(.6)
#    L_diff = L_expected - L_estimated
#
#    for i in range(0,L_diff.shape[0]):
#        for j in range(0,L_diff.shape[1]):
#            pass
#
#
#    fh.test_linking(A,L_expected)
#    pl.matshow(A)
#    pl.title("Assignment")
#    pl.matshow(L_expected)
#    pl.title("Expected")
#    pl.matshow(L_estimated)
#    pl.title("Estimated")
#    pl.matshow(L_diff)
#    pl.title("Diff")
#    pl.colorbar()
#    pl.show()

    
    
    
