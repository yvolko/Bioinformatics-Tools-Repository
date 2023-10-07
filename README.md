# Bioinformatics-Tools-Repository
Collection of bioinformatic tools to work with fastq sequences, DNA and RNA sequences and protein sequences
# Protein Info

This tool supports standard 20 amino acids. Any modifications of amino acids are not supported. You can write amino acids in any case (lower, upper or mixed). 
This project consists of one function "protein_analysis" that helps user to:
- predict molecular weight of amino acid (aa) sequences
- translate aa sequences from one-letter to three-letter code
- calculate total amount of each amino acid in the sequences
- make DNA based codon optimization for the introduced amino acid sequences with the support for 3 cell types: Esherichia coli, Pichia pastoris, Mouse
- calculate length of amino acid sequences
- count the number of atoms of each type in a sequence (brutto formula)  <br/>

Tool is coded with Python.

## How to use:
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

## List of procedures:

- `molecular_weight` — returns list of float values, that indicate predicted molecular weights of given aa sequences (in kDa)
- `one_letter_to_three` — will return list of strings, containing the same sequences written in three-letter code
- `get_amino_acid_sum` — сounts the amount of each amino acid in the injected protein sequences
- `codon_optimization` — makes codon-optimized DNA based on the introduced amino acid sequences for 3 types of cells: Esherichia coli, Pichia pastoris, Mouse
- `length` — calculates length of amino acid sequences 
- `brutto_count` — counts the number of atoms of each type in a sequence

## Example of use:

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

## Private policy and contacts
This tool can be freely distributed and used.
<br/>
If you have any suggestions for improving the tool or if you find a bug, please contact us by email.
<br/>
This tool was developed by the "workaholics" team:
<br/>
Yulia Volkova volkova.yulia.leonidovna@gmail.com
<br/>
Dasha Sokolova kalabanova_dasha@mail.ru
<br/>
Team leader: Ivan Kozin ivan.d.kozin@gmail.com
<br/>
Team photo:
![Снимок экрана 2023-09-29 210559_2](https://github.com/ivandkoz/HW4_Functions2_Kozin/assets/63678919/ad1302a1-d139-4c82-b7eb-d5b9ac1897e8)

## Personal contribution
`Ivan Kozin` (team leader) worte functions:
- length
- brutto_count
- is_amino_acid
- name_transform
- is_length_divisible_by_3
- is_amino_acid_three_letter
- managed work with guthub repository

`Dasha Sokolova` (co-leader) wrote functions: 
- get_amino_acid_sum
- codon_optimization functions
  
`Yulia Volkova` (co-leader) wrote functions:
- main (protein_analysis)
- molecular_weight
- one_letter_to_three functions
  
Writting README, debugging code and testing it has been done by the efforts of all team.



