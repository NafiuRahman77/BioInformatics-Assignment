import sys
import csv
file = sys.argv[1]

# def score(motifs, k):
#     consensus = ''
#     score = 0
#     for i in range(k):
#         col = [motif[i] for motif in motifs]
#         max_freq = max(col, key=col.count)
#         consensus += max_freq
#         score += sum(1 for base in col if base != max_freq)
#     return score, consensus

def score_h (consensus, motif, k):
    score = 0
    for i in range(k):
        if consensus[i] != motif[i]:
            score += 1
    return score

def readfile():
    motifs = []
    scores = []
    file_path = file
    target_motif_id = '1'
    consensus = ''

    with open(file_path, 'r', encoding='utf-8') as tsv_file:
        tsv_reader = csv.DictReader(tsv_file, delimiter='\t')
        first = False
        for row in tsv_reader:    
            if first == False:
                consensus = row['motif_ID'].split('-')[1]
            first = True
            if row['motif_ID'][0] == '#':
                break      
            if row['motif_ID'][0] == target_motif_id:
                site_seq = row['site_Sequence']
                idx = int(row['seq_ID'])
                if idx > len(motifs):
                    motifs.extend([''] * (idx + 1 - len(motifs)))
                    scores.extend([0] * (idx + 1 - len(scores)))
                new_score = score_h(consensus, site_seq, len(consensus))
                if (scores[idx-1] == 0 or scores[idx-1] > new_score):
                    scores[idx-1] = new_score
                    motifs[idx-1] = row['site_Sequence']
        print (motifs)
        print (scores)
    return motifs, len(motifs[0]), consensus, sum(scores)

def main():
    motifs,k, consensus, score = readfile()
    print(motifs)
    print(consensus)
    print(score)

if __name__ == '__main__':

    main()

