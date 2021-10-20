"""Milestone1 Script
"""
import itertools
from functools import reduce

genetic_code = {
    'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L', 'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
    'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M', 'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
    'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S', 'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'UAU': 'Y', 'UAC': 'Y', 'UAA': 'Stop', 'UAG': 'Stop', 'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'UGU': 'C', 'UGC': 'C', 'UGA': 'Stop', 'UGG': 'W', 'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}
 

#Haocen Xie
def s(dna):
    Dict = {}
    count_A = count_C = count_G = count_T = 0
    dna = dna.upper()
    for i in dna:
        if i == 'A':
            count_A += 1
            Dict['A'] = count_A
        elif i == 'C':
            count_C += 1
            Dict['C'] = count_C
        elif i == 'G':
            count_G += 1
            Dict['G'] = count_G
        elif i == 'T':
            count_T += 1
            Dict['T'] = count_T
        
    return Dict


#Haocen Xie
def dna2rna(dna):
    rna = ''
    for letter in dna:
        if letter == 'T':
            rna += 'U'
        else:
            rna += letter
    return rna


#Haocen Xie
def reverse_complement(dna):
    reverse = ''
    for i in dna:
        if i == 'A':
            i = 'T'
        elif i == 'T':
            i = 'A'
        elif i == 'C':
            i = 'G'
        elif i == 'G':
            i = 'C'
        reverse = i + reverse
    return reverse


#Haocen Xie
def mendels_law(hom,het,rec):
    sum = hom + het + rec
    reces = (rec/sum)*((rec-1)/(sum-1))
    heter = (het/sum)*(het-1)/(sum-1)
    heter_reces = (rec/sum)*(het/(sum-1))+(het/sum)*(rec/(sum-1))
    result = 1 - (reces + heter * 1/4+ heter_reces * 1/2 )
    return result


#Haocen Xie
def fibonacci_rabbits(n, k):
    f_1, f_2 = 1, 1
    for i in range(n - 1):
        f_2, f_1 = f_1, f_1 + (f_2 * k)
    return f_2


#Manxin Chen
def GC_content(dna_list):
    stats = [
    (seq.count('C')+seq.count('G'))/len(seq)                    
    for seq in dna_list                                         
    ]
    index = max(range(len(stats)), key = lambda i:stats[i])     
    return (index, 100*stats[index])


#Manxin Chen
def rna2codon(rna):
    codon = []                                                  
    for i in range(0, len(rna), 3):                             
    	if genetic_code[rna[i:i+3]] == 'Stop':                  
    		break
    	codon.append(genetic_code[rna[i:i+3]])                  
    return ''.join(codon)


#Manxin Chen
def locate_substring(dna_snippet, dna):                     
    return [
        i for i in range(len(dna)-len(dna_snippet))             
        if dna[i:i+len(dna_snippet)] == dna_snippet             
    ]



#Manxin Chen
def hamming_dist(dna1, dna2):
    return sum( dna2[i] != dna1[i] for i in range(0,len(dna1))) 


#Manxin Chen
def source_rna(protein):
    rna = []                                                    
    for i in protein:                                          
        rna.append(
            [v for v,k in genetic_code.items()                 
            if k==i]                                            
            )
    rna.append(
        [v for v,k in genetic_code.items() if k == 'Stop']      
        )
    return reduce(lambda x,y: x*y,[len(i) for i in rna])        


#Ali
def count_dom_phenotype(genotypes):
    offspring = 0

    for i in range(len(genotypes)):

        if i == 0:
            offspring = (genotypes[0]*2*1)
        elif i == 1:
            offspring += (genotypes[1]*2*1)
        elif i == 2:
            offspring += (genotypes[2]*2*1)
        elif i == 3:
            offspring += (genotypes[3]*2*0.75)
        elif i == 4:
            offspring += (genotypes[4]*2*0.5)
        elif i == 5:
            offspring += (genotypes[4]*2*0)
        else:
            return print('String is too long')
    return offspring


#Ali
def splice_rna(dna, intron_list):
    for i in intron_list:
        dna = dna.replace(i,"")
        
    rna = dna.replace("T","U")
    return rna2codon(rna)