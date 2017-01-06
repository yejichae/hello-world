
from __future__ import division
from collections import defaultdict

input_file = 'practice6.fasta'

def import_file(input_file):
    with open(input_file) as f_in:
        txt = (line.rstrip() for line in f_in)
        txt = list(line for line in txt if line)
    return txt

seq = import_file(input_file)

D_fasta = defaultdict(str)
for line in seq:
    if line.startswith('>'):
        scaffold_name = line.replace('>', '')
    else:
        D_fasta[scaffold_name] += line

def genome_seq(seq):
    genome_seq = ''
    for seq1 in D_fasta.values():
        genome_seq += seq1
    return genome_seq

genome_seq = genome_seq(seq)

def gc_content(seq):
    gc_count = seq.count('G') + seq.count('C')
    gc_content = gc_count / len(seq) * 100
    return str(gc_content)

print '''1. GC content of whole genome: ''' + gc_content(genome_seq)
print '''2. Genome size(bp): ''' + str(len(genome_seq))

def longest_scaffold(seq):
    longest = ''
    longest_scaffold = ''
    for scaffold_name in D_fasta.keys():
        if len(D_fasta[scaffold_name]) > len(longest):
            longest = D_fasta[scaffold_name]
            longest_scaffold = scaffold_name
    return longest_scaffold
          
print '''3. Longest scaffold(name&length): ''' + longest_scaffold(seq) + ' & ' + str(len(D_fasta[longest_scaffold(seq)]))
         
def average_scaffold_size(seq):
    size_average = len(genome_seq) / len(D_fasta.values())
    return size_average     

print '''4. Average scaffold size: ''' + str(average_scaffold_size(seq))
       
def scaffold_size(seq):
    scaffold_size = []
    for scaffold in D_fasta.values():
        scaffold_size.append(len(scaffold))
    scaffold_size.sort(reverse=True)
    return scaffold_size
    
def N50(seq):
    half_sum = len(genome_seq) / 2
    making_half = 0
    for i in range(len(scaffold_size(seq))):
        if half_sum > making_half:
            making_half += scaffold_size(seq)[i]
        else:
            break
    return scaffold_size(seq)[i-1]

print '''5. N50: ''' + str(N50(seq))

def number_more_1kb(seq):
    count = 0
    for size in scaffold_size(seq):
        if size > 1000:
            count += 1
    return count

print '''6. The number of scaffolds of which length is more than 1 kb: ''' + str(number_more_1kb(seq))

print '''7. Number of scaffolds: ''' + str(len(scaffold_size(seq)))
