import os
from dataclasses import dataclass


@dataclass
class FastaRecord:
    fasta_id: str
    description: str
    seq: str

    def __repr__(self):
        return f"FastaRecord:\nid = {self.fasta_id} \n" \
               f"description = {self.description} \n" \
               f"sequence = {self.seq}"


class OpenFasta:
    def __init__(self, file_path, mode='r'):
        self.file_path = file_path
        self.mode = mode
        self.file = None
        self._current_line = None
        self._id = None
        self._description = None

    def __enter__(self):
        self.file = open(self.file_path, self.mode)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.file.close()

    def __iter__(self):
        return self

    def __next__(self):
        record = self.read_record()
        if record is None:
            raise StopIteration
        return record

    def read_record(self):
        if self._current_line is None:
            self._current_line = self.file.readline().strip()
        if self._current_line.startswith(">"):
            self._id, self._description = self._current_line.split(" ", 1)
        if not self._current_line:
            return None

        seq = ""
        fasta_record = FastaRecord(fasta_id=self._id[1:],
                                   description=self._description,
                                   seq="")
        while True:
            self._current_line = self.file.readline().strip()
            if not self._current_line:
                break
            if self._current_line.startswith(">"):
                self._id, self._description = self._current_line.split(" ", 1)
                break
            seq += self._current_line
        fasta_record.seq = seq
        return fasta_record

    def read_records(self):
        records = []
        while True:
            record = self.read_record()
            if record:
                records.append(record)
            else:
                break
        return records


def convert_multiline_fasta_to_oneline(input_fasta: str,
                                       output_fasta: str = '') -> None:
    """
    Converts fasta file where sequences are written on multiple lines
    to fasta file where each sequence is written on one line
    :param input_fasta: string, valid path to the input file
    with the name of the file and file extension at the end
    :param output_fasta: string, name of the output file.
    By default - input fasta file name is taken and "_oneline" is added to it.
    :return: None. Will write output to the file.
    """
    if output_fasta == '':
        output_fasta = os.path.basename(input_fasta)
        output_fasta = output_fasta.split('.')
        output_fasta = f'{output_fasta[0]}_oneline.{output_fasta[1]}'

    first_line = True
    with open(input_fasta, 'r', encoding='utf-8') as input_file:
        try:
            with open(output_fasta, 'x', encoding='utf-8') as output_file:
                for line in input_file.readlines():
                    if line.startswith('>'):
                        if first_line:
                            output_file.write(line)
                            first_line = False
                        else:
                            output_file.write('\n')
                            output_file.write(line)
                    else:
                        output_file.write(line.strip())
        except FileExistsError:
            print(f'File with the {output_fasta} name already exist.')
            print('Please use another name.')


def select_genes_from_gbk_to_fasta(input_gbk: str,
                                   genes: str,
                                   n_before: int = 1,
                                   n_after: int = 1,
                                   output_fasta: str = '') -> None:
    """
    Select neighbours of the gene(s) provided from gbk input file.
    :param input_gbk: string, valid path to the input gbk file
    with the name of the file and file extension at the end
    :param genes: string, genes of interest separated with space
    :param n_before: integer, how many genes to take before gene of interest
    :param n_after: integer, how many genes to take after gene of interest
    :param output_fasta: string, name of the output fasta file.
    By default - input gbk file name is taken and "_selected_genes"
    is added to it.
    :return: None. Will write output to the file.
    """
    if output_fasta == '':
        output_fasta = os.path.basename(input_gbk)
        output_fasta = output_fasta.split('.')
        output_fasta = f'{output_fasta[0]}_selected_genes'

    genes = genes.split()
    input_dict_position = {}
    input_dict_sequences = {}

    with open(input_gbk, 'r', encoding='utf-8') as input_file:
        gene_order_number = 1
        line = input_file.readline()
        while line:
            if 'CDS' in line:
                line = input_file.readline()
                if '/gene' in line or '/locus_tag' in line:
                    name = line.split('\"')[1]
                    input_dict_position[name] = gene_order_number
                    gene_order_number += 1
                    while '/translation' not in line:
                        line = input_file.readline()
                    input_dict_sequences[name] = line.split('\"')[1].strip()
                    if line.strip()[-1] == '\"':
                        pass
                    else:
                        line = input_file.readline().strip()
                        line = line.split('\"')[0].strip()
                        input_dict_sequences[name] += line
                        while '\"' not in line:
                            line = input_file.readline()
                            line.split('\"')[0].strip()
                            input_dict_sequences[name] += line
            line = input_file.readline()

    neighbours_of_gene = {}
    for gene in genes:
        if gene in input_dict_position.keys():
            current_gene_position = input_dict_position[gene]
            left_border = current_gene_position - n_before
            right_border = current_gene_position + n_after
            left_border = max(left_border, 0)
            if right_border >= len(input_dict_position):
                right_border = len(input_dict_position)
            for key, value in input_dict_position.items():
                if value in range(left_border,
                                  right_border + 1) and key != gene:
                    neighbours_of_gene[key] = input_dict_sequences[key]
        else:
            raise ValueError(f'The gene {gene} is not in the data.')

    try:
        with open(f'{output_fasta}.fasta', mode='x', encoding='utf-8') as file:
            for key, value in neighbours_of_gene.items():
                file.write(f'>{key}\n')
                file.write(f'{value}\n')
    except FileExistsError:
        print('File with the provided name already exist.')
        print('Please use another name.')


def change_fasta_start_pos(input_fasta: str,
                           shift: int,
                           output_fasta: str) -> None:
    """
    Will return fasta read from (shifted to) the position equals to shift.
    :param input_fasta: string, valid path to the input file
    (this file should contain only one sequence:
    first line starts with '>', second line - sequence)
    :param shift: integer, step to shift start of reading the sequence
    :param output_fasta: string, name of the output file
    :return: None, will print data to fasta file
    """
    output_fasta = output_fasta + '.fasta'
    with open(input_fasta, 'r', encoding='utf-8') as input_file:
        line_1 = input_file.readline()
        line_2 = input_file.readline().strip()
        shifted_line = line_2[shift:] + line_2[:shift]
        try:
            with open(output_fasta, 'x', encoding='utf-8') as output_file:
                output_file.write(line_1)
                output_file.write(shifted_line)
        except FileExistsError:
            print('File with the provided name already exist.')
            print('Please use another name.')
