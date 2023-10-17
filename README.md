# Bioinformatics-Tools-Repository
Collection of bioinformatic tools to work with fastq format, gbk, DNA and RNA sequences and protein sequences. 
# Working with DNA and RNA sequences

This tool helps to get reverse, complement, reverse-complement and transcribed sequences. As well as define if sequense is DNA or is RNA. 

# Working with protein sequences

This tool supports standard 20 amino acids. Any modifications of amino acids are not supported. You can write amino acids in any case (lower, upper or mixed). 
1. Alanine (A, Ala)
2. Arginine (R, Arg)
3. Asparagine (N, Asn)
4. Aspartic Acid (D, Asp)
5. Cysteine (C, Cys)
6. Glutamine (Q, Gln)
7. Glutamic Acid (E, Glu)
8. Glycine (G, Gly)
9. Histidine (H, His)
10. Isoleucine (I, Ile)
11. Leucine (L, Leu)
12. Lysine (K, Lys)
13. Methionine (M, Met)
14. Phenylalanine (F, Phe)
15. Proline (P, Pro)
16. Serine (S, Ser)
17. Threonine (T, Thr)
18. Tryptophan (W, Trp)
19. Tyrosine (Y, Tyr)
20. Valine (V, Val)
    
This project consists of one function "protein_analysis" that helps user to:
- predict molecular weight of amino acid (aa) sequences
- translate aa sequences from one-letter to three-letter code
- calculate total amount of each amino acid in the sequences
- make DNA based codon optimization for the introduced amino acid sequences with the support for 3 cell types: Esherichia coli, Pichia pastoris, Mouse
- calculate length of amino acid sequences
- count the number of atoms of each type in a sequence (brutto formula)  <br/>

# Working with fastq format
This tools choses sequences that satisfy conditions of sertain length, GC content and sequncing quality (phred33 standard). <br/>
Learn more about fastq files [here](https://stepik.org/lesson/32398/step/1?unit=12379). <br/>
Learn more about quality score encoding in fastq files [here](https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/QualityScoreEncoding_swBS.htm).

# Convert multiline fasta file to oneline fasta file
Converts fasta file where sequences are written on multiple lines to fasta file where each sequence is written as single line.

# Select genes from gbk to fasta
Selects neighbouring genes to the genes of interest.

## How to use:
### DNA/RNA tool
**run_dna_rna_tools**(**args, procedure) <br/>
**Parametrs:**
> ***args** : **sequence of str** <br/>
> &nbsp;&nbsp;&nbsp;&nbsp;Any number of lines with DNA only or RNA only sequences <br/>
  **procedure** : ***str*** <br/>
> &nbsp;&nbsp;&nbsp;&nbsp;The name of the operation you want to perform. The following types of procedures are supported: <br/>
>>  
>> - ***transcribe***: transcribes DNA sequnces
>> - ***reverse***: return reverse DNA or RNA sequence
>> - ***complement***: return complement DNA sequence
>> - ***reverse_complement***: return reverse-complement DNA sequence
>> - ***get_nucl_acid_type***: checks wether give sequence is DNA or RNA (ND-not defined)

Call the "run_dna_rna_tools" funcion with following arguments.
Requred arguments:
- sequence of strings representing only DNA or only RNA. Please do not use a mixture of DNA/RNA sequences in the same function call!
- name of procedure as string (see list of precedures)

### Protein tool
**protein_analysis**(**args, procedure, cell_type=None, letter_format=1*) <br/>
**Parametrs:**
> ***args** : **sequence of str** <br/>
> &nbsp;&nbsp;&nbsp;&nbsp;Any number of lines with amino acid sequences <br/>
    **procedure** : ***str*** <br/>
> &nbsp;&nbsp;&nbsp;&nbsp;The name of the operation you want to perform. The following types of procedures are supported: <br/>
>>  
>> - ***molecular_weight***: calculates predicted molecular weight of amino acid sequences in kDa
>> - ***one_letter_to_three***: translate aa sequences from one-letter to three-letter code
>> - ***get_amino_acid_sum***: calculates total amount of each amino acid in the sequences
>> - ***codon_optimization***: makes DNA based codon optimization for the introduced amino acid sequences, support 3 types of cells. Can only be used in conjunction with **cell_type**: `Esherichia coli`, `Pichia pastoris`, `Mouse`
>> - ***length***: calculates length of amino acid sequences 
>> - ***brutto_count***: counts the number of atoms of each type in a sequence
>> 
>    **cell_type** : ***str, defalut None*** <br/>
> &nbsp;&nbsp;&nbsp;&nbsp;The type of cells for which optimization is applied. Cell types supported:<br/>
>>
>> - `Esherichia coli` *or* `E.coli`
>> - `Pichia pastoris` *or* `P.pastoris`
>> - `Mouse` *or* `mouse`
>> 
>    **letter_format** : ***int, defalut 1*** <br/>
> &nbsp;&nbsp;&nbsp;&nbsp;Specifies the format for receiving amino acid sequences. Either one-letter (**letter_format** = 1) or three-letter sequences (**letter_format** = 3) <br/>
>

Call the "protein_analysis" funcion with following arguments.
Requred arguments:
- tuple of protein sequences written one letter or three letter code without stop codos. Please do not use sequences in different formats in the same function call!
- name of procedure as string (see list of precedures)
- format of code for the protein sequences as int: 1 for one letter, 3 for three letter code
Optional argument:
- cell type (required only for codon_optimization procedure). Accepted cell types Esherichia coli, Pichia pastoris, Mouse

### FASTQ tool
**fastq_thresholding**(input_path, output_filename = '', gc_bounds = (0, 100), length_bounds = (0, 2 ** 32), quality_threshold = 0) <br/>
**Parametrs:**
> ***input_path** : **str** <br/>
> &nbsp;&nbsp;&nbsp;&nbsp;String, valid path to the input file with the name of file and extension of the file <br/> <br/>
> ***output_filename** : **str** <br/>
> &nbsp;&nbsp;&nbsp;&nbsp;String, name of the output file. By default - input fasta file name is taken. <br/> <br/>
    **gc_bounds** : ***tuple, defalut (0, 100)*** <br/>
> &nbsp;&nbsp;&nbsp;&nbsp;Tuple: first value is lower boundary of GC-content and second value - upper boundry. Function will filter sequences that are in between given values (including given values). If only one number is given it is taken as upper bound and lower bound is taken as 0. <br/> <br/>
    **length_bounds** : ***tuple, defalut (0, 2 power(32))*** <br/>
> &nbsp;&nbsp;&nbsp;&nbsp;Tuple: first value is lower boundary of length and second value - upper boundry. Function will filter sequences that are in between given values (including given values). If only one number is given it is taken as upper bound and lower bound is taken as 0. <br/> <br/>
    **quality_threshold** : ***int, defalut 1*** <br/>
> &nbsp;&nbsp;&nbsp;&nbsp;If average sequnce quality of given sequence is above or equal to given number sequence will be selected. <br/>
>

Call the "fastq_thresholding" funcion with following arguments. <br/>

Requred arguments: <br/>
- input_path <br/> 

Optional argument: <br/>
- output_filename
- gc_bounds
- length_bounds
- quality_threshold

### Convert multiline fasta file to oneline fasta file
**convert_multiline_fasta_to_oneline**(input_fasta, output_fasta = '') <br/>
**Parametrs:**
> ***input_fasta** : **str** <br/>
> &nbsp;&nbsp;&nbsp;&nbsp;String, valid path to the input file with the name of file and extension of the file <br/> <br/>
> ***output_fasta** : **str** <br/>
> &nbsp;&nbsp;&nbsp;&nbsp;String, name of the output file. By default - input fasta file name is taken and "_oneline" is added to it. <br/> <br/>

Call the "convert_multiline_fasta_to_oneline" funcion with following arguments. <br/>

Requred arguments: <br/>
- input_fasta <br/> 

Optional argument: <br/>
- output_fasta

### Select genes from gbk to fasta
**select_genes_from_gbk_to_fasta**(input_gbk, genes, n_before = 1, n_after = 1, output_fasta = '') <br/>
**Parametrs:**
> ***input_gbk** : **str** <br/>
> &nbsp;&nbsp;&nbsp;&nbsp;String, valid path to the input file with the name of file and extension of the file <br/> <br/>
> ***genes** : **str** <br/>
> &nbsp;&nbsp;&nbsp;&nbsp;String, genes of interest separated with space. <br/> <br/>
> ***n_before** : **int** <br/>
> &nbsp;&nbsp;&nbsp;&nbsp;Integer, how many genes to select before gene of interest. By default equal to 1. <br/> <br/>
> ***n_after** : **int** <br/>
> &nbsp;&nbsp;&nbsp;&nbsp;Integer, how many genes to select after gene of interest. By default equal to 1. <br/> <br/>
> ***output_fasta** : **str** <br/>
> &nbsp;&nbsp;&nbsp;&nbsp;String, name of the output fasta file. By default name of the input gbk file is taken and "_selected_genes" is added. <br/> <br/>

Call the "select_genes_from_gbk_to_fasta" funcion with following arguments. <br/>

Requred arguments: <br/>
- input_gbk <br/>
- genes <br/>

Optional argument: <br/>
- n_before
- n_after
- output_fasta

## List of procedures:
- `transcribe` — returns list of strings with transcribed sequences (or one string if one sequence is given)
- `reverse` — returns list of strings with reversed sequences (or one string if one sequence is given)
- `complement` — returns list of strings with complement sequences (or one string if one sequence is given)
- `reverse_complement` — returns list of strings with reverse-complement sequences (or one string if one sequence is given)
- `get_nucl_acid_type` — returns list of string with `DNA`, `RNA` or `ND` ()not defined in it (or one string if one sequence is given)
- `molecular_weight` — returns list of float values, that indicate predicted molecular weights of given aa sequences (in kDa)
- `one_letter_to_three` — will return list of strings, containing the same sequences written in three-letter code
- `get_amino_acid_sum` — сounts the amount of each amino acid in the injected protein sequences
- `codon_optimization` — makes codon-optimized DNA based on the introduced amino acid sequences for 3 types of cells: Esherichia coli, Pichia pastoris, Mouse
- `length` — calculates length of amino acid sequences 
- `brutto_count` — counts the number of atoms of each type in a sequence

## Example of use:

```python
run_dna_rna_tools('ATG', 'transcribe') # 'AUG'
run_dna_rna_tools('ATG', 'reverse') # 'GTA'
run_dna_rna_tools('AtG', 'complement') # 'TaC'
run_dna_rna_tools('ATg', 'reverse_complement') # 'cAT'
run_dna_rna_tools('ATG', 'aT', 'reverse') # ['GTA', 'Ta']
run_dna_rna_tools('AUG', 'aT', 'aCg', 'is_dna_is_rna') # ['RNA', 'DNA', 'ND']
```

```python
protein_analysis("ACD", "AD", procedure="one_letter_to_three", letter_format=1) # ['AlaCysAsp', 'AlaAsp']
protein_analysis("AlaAspLys", "AlaAsp", procedure="molecular_weight", letter_format=3) # [0.37, 0.22]
protein_analysis("ACD", "AD", procedure="get_amino_acid_sum") # [{'A': 1, 'C': 1, 'D': 1, 'E': 0, 'F': 0, 'G': 0, 'H': 0, 'I': 0, 'K': 0, 'L': 0, 'M': 0, 'N': 0, 'P': 0, 'Q': 0, 'R': 0, 'S': 0, 'T': 0, 'V': 0, 'W': 0, 'Y': 0},
                                                                        # {'A': 1, 'C': 0, 'D': 1, 'E': 0, 'F': 0, 'G': 0, 'H': 0, 'I': 0, 'K': 0, 'L': 0, 'M': 0, 'N': 0, 'P': 0, 'Q': 0, 'R': 0, 'S': 0, 'T': 0, 'V': 0, 'W': 0, 'Y': 0}]
protein_analysis("ACD", "AD", procedure="codon_optimization", cell_type = 'E.coli', letter_format=1) # ['GCGTGCGAT', 'GCGGAT']
protein_analysis("acDEFGHIKLMNPQRSTVwy", "ad", procedure="length", letter_format=1) # [20, 2]
protein_analysis("FGHIKLMNPQ", "PQRSTVwy", "adN", procedure="brutto_count", letter_format=1)
# [{'C': 54, 'H': 103, 'N': 15, 'O': 22, 'S': 1}, {'C': 48, 'H': 83, 'N': 23, 'O': 18, 'S': 3}, {'C': 11, 'H': 22, 'N': 4, 'O': 9, 'S': 0}]
```


```python
fastq_thresholding("C:\\Users\\name\\Desktop\\pythonProject\\example_data.txt", "my_file22", length_bounds=(1,1000)) # in windows OS
fastq_thresholding("/Users/volko/Desktop/pythonProject/example_data.txt", "my_file22", length_bounds=(1,1000)) # in Linux/Mac OS
```

```python
select_genes_from_gbk_to_fasta('C:\\Users\\name\\Desktop\\pythonProject\\example_gbk.gbk', "rrrD_1 pxpB", n_before=2, n_after=2) # in windows OS
fastq_thresholding("/Users/name/Desktop/pythonProject/example_gbk.gbk", "my_file22", length_bounds=(1,1000)) # in Linux/Mac OS
```

```python
convert_multiline_fasta_to_oneline('C:\\Users\\name\\Desktop\\pythonProject\\example_multiline_fasta.fasta')) # in windows OS
convert_multiline_fasta_to_oneline("/Users/name/Desktop/pythonProject/example_multiline_fasta.fasta", "my_file22", length_bounds=(1,1000)) # in Linux/Mac OS
```

## Input requirements and possible errors:
 - **It is important to indicate the type of operation. An error occurs when you enter an incorrect operation type**
```python
protein_analysis("FGHIKLMNPQ", "PQRSTVwy", "adN", procedure="brutto", letter_format=1)
# ValueError: Requested procedure is not defined
```
 - **To perform the coden_optimization operation, you must enter cell_type (None by default). Otherwise an error message is displayed**
```python
protein_analysis('AlaCysAsp', 'AlaAsp', procedure="codon_optimization", cell_type='Rat', letter_format=3) 
# ValueError: Type Rat is not supported. The following types of organisms are available for codon optimization: Esherichia coli, Pichia pastoris, Mouse
```
 - **By default, entering amino acid sequences in a single-letter format in any case is supported. To enter in three-letter format in any case, you need to specify letter_format = 3. <br/> If an unknown format is entered, an error message is displayed.**
```python
protein_analysis("ACD", "AD", procedure="one_letter_to_three", cell_type='E.coli', letter_format=2)
# ValueError: Error unsupported letter_format. Only letter_formats 1 and 3 are supported
```
 - **If letter_format = 1 is specified, but all sequences are similar to the three-letter amino slot encoding, a notification will be displayed warning**
```python
protein_analysis("LYSlys", "HishisHis", procedure="get_amino_acid_sum", letter_format=1)
# Warning: all your sequences are similar to three-letter ones. Check the letter_format value
```
 - **If a single-letter amino acid input format is specified, but at least one amino acid slot is not standard or is written incorrectly, an error message is displayed**
```python
protein_analysis("BBB", procedure="get_amino_acid_sum", letter_format=1))
# ValueError: Error B is not an amino acid. Correct your input
```
 - **If a three-letter amino acid input format is specified, but at least one amino acid slot is not standard or is written incorrectly, an error message is displayed**
```python
protein_analysis("Al", procedure="get_amino_acid_sum", letter_format=3)
# ValueError: Error al is incorrect form of amino acid notation. Correct your input
protein_analysis("AluLysArg", procedure="get_amino_acid_sum", letter_format=3)
# ValueError: Error alu is not an amino acid. Correct your input
```

 - **If file with the name equal to output_filename (in fastq_thresholding function) or equal output_fasta (in convert_multiline_fasta_to_oneline or in select_genes_from_gbk_to_fasta functions) already exists FileExistsError will occure**
```python
fastq_thresholding("C:\\Users\\name\\Desktop\\pythonProject\\example_data.txt", "my_file22", length_bounds=(1,1000))
# FileExistsError: File with the provided name already exist. Please use another name.
convert_multiline_fasta_to_oneline("C:\\Users\\name\\Desktop\\pythonProject\\example_multiline_fasta.fasta", "my_file22", )
# FileExistsError: File with the provided name already exist. Please use another name.
select_genes_from_gbk_to_fasta("C:\\Users\\name\\Desktop\\pythonProject\\example_gbk.gbk", "my_file22", "rrrD_1 pxpB", n_before=2, n_after=2)
# FileExistsError: File with the provided name already exist. Please use another name.
```

 - **If at least one of the genes in the parameter genes (function select_genes_from_gbk_to_fasta) is not in the input gbk file ValueError will occure**
```python
select_genes_from_gbk_to_fasta('C:\\Users\\name\\Desktop\\pythonProject\\example_gbk.gbk', "AMMMttt", n_before=2, n_after=2)
# ValueError: 'The gene AMMMttt is not in the data.'
```

## Personal contribution
Protein tool was written in team with:
`Ivan Kozin` wrote functions:
- length
- brutto_count
- is_amino_acid
- name_transform
- is_length_divisible_by_3
- is_amino_acid_three_letter
- managed work with guthub repository

`Dasha Sokolova` wrote functions: 
- get_amino_acid_sum
- codon_optimization functions
  
Everyting else has been written by `Yulia Volkova`.

All tools are coded with Python.



