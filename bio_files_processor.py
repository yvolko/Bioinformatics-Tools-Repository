import os
import codecs

def convert_multiline_fasta_to_oneline(input_fasta: str, output_fasta: str = ''):
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
