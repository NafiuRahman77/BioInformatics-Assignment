import csv
def score(motifs, k):
    consensus = ''
    score = 0
    for i in range(k):
        col = [motif[i] for motif in motifs]
        max_freq = max(col, key=col.count)
        consensus += max_freq
        score += sum(1 for base in col if base != max_freq)
    return score, consensus

def readfile():
    motifs = []
    file_path = 'streme/sites-3.tsv'
    target_motif_id = '1'
    
    with open(file_path, 'r', encoding='utf-8') as tsv_file:
        tsv_reader = csv.DictReader(tsv_file, delimiter='\t')
        for row in tsv_reader:          
            if row['motif_ID'][0] == target_motif_id:
                #print(row['site_Sequence'])
                motifs.append(row['site_Sequence'])

    return motifs, len(motifs[0])

def main():
    motifs,k = readfile()
    print(motifs)
    print(score(motifs, k))

if __name__ == '__main__':

    main()

