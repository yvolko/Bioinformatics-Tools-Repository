import os
import codecs

def convert_multiline_fasta_to_oneline(input_fasta: str, output_fasta: str = '') -> None:
    if output_fasta == '':
        output_fasta = os.path.basename(input_fasta)
        output_fasta = output_fasta.split('.')
        output_fasta = f'{output_fasta[0]}_oneline.{output_fasta[1]}'

    with open (input_fasta) as input_file:
        try:
            with open(output_fasta, 'x') as output_file:
                output_file.write(input_file.readline())
                for line in input_file.readlines():
                    if line.startswith('>'):
                        output_file.write('\n')
                        output_file.write(line)
                    else:
                        output_file.write(line.strip())
        except FileExistsError:
            print('File with the provided name already exist. Please use another name.')


def select_genes_from_gbk_to_fasta(input_gbk: str, genes: str, n_before: int = 1,
                                   n_after: int = 1, output_fasta: str = '') -> None:
    if output_fasta == '':
        output_fasta = os.path.basename(input_gbk)
        output_fasta = output_fasta.split('.')
        output_fasta = f'{output_fasta[0]}_selected_genes'

    genes = genes.split()
    input_dict_position = {}
    input_dict_sequences = {}

    with codecs.open(input_gbk, 'r', 'utf-8') as input_file:
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
                        break
                    else:
                        line = input_file.readline().strip()
                        while 'CDS' not in line:
                            new_value = input_dict_sequences[name] + line.strip()
                            input_dict_sequences[name] = new_value
                            line = input_file.readline().replace('\"', '')
            line = input_file.readline()

    for gene in genes:
        if gene in input_dict_position.keys():
            current_gene_position = input_dict_position[gene]
            neighbours_of_gene = {}
            left_border = current_gene_position - n_before
            right_border = current_gene_position + n_after
            if left_border < 0:
                left_border = 0
            if right_border >= len(input_dict_position):
                right_border = len(input_dict_position)
            for key, value in input_dict_position.items():
                if value in range (left_border, right_border+1) and key != gene:
                    neighbours_of_gene[key] = input_dict_sequences[key]
        else:
            raise ValueError(f'The gene {gene} is not in the data.')

    try:
        with open(f'{output_fasta}.fasta', mode='x') as file:
            for key, value in neighbours_of_gene.items():
                file.write(f'>{key}\n')
                file.write(value + '\n')
    except FileExistsError:
        print('File with the provided name already exist. Please use another name.')
