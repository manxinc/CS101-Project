"""Milestone2 Script
"""
# Chen
def find_splice(dna_motif,dna):
    j = 0
    splice = []
    for i in range(len(dna)):
        if j >= len(dna_motif):
            break
        if dna[i] == dna_motif[j]:
            splice.append(i)
            j +=1
    if [dna[j] for j in splice] != [i for i in dna_motif]:
        splice = []
    return splice


# Chen
def shared_motif(dna_list):
    
    s = dna_list[0]
    shared = ''
    
    for i in range(len(s)):
        for j in range(i+1,len(s)+1):
            stem = s[i:j]
            k = 1
            for k in range(1,len(dna_list)):
                if stem not in dna_list[k]:
                    break
                if (k+1 == len(dna_list) and len(shared) < len(stem)):
                    shared = stem
                
    return shared

#Chen
def get_edges(DNA_dict):
    edges = []
    dict_reverse = {v:k for k,v in DNA_dict.items()}
    for i in DNA_dict.values():
        for j in DNA_dict.values():
            if i == j:
                continue
            if i[-3:] == j[0:3]:
                edges.append((dict_reverse[i],dict_reverse[j]))
    return edges