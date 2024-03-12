import random

def read_sequences(filepath):
    with open(filepath, 'r') as file:
        sequences = [line.strip() for line in file.readlines()]
    return sequences

def randomly_select_kmers(sequences, k):
    return [seq[start:start+k] for seq in sequences for start in [random.randint(0, len(seq) - k)]]

def create_profile(motifs, k):
    profile_matrix = [[1 for _ in range(k)] for _ in range(4)]
    nucleotides = 'ACGT'
    for index, motif in enumerate(motifs):      
        for position, nucleotide in enumerate(motif):
            row = nucleotides.index(nucleotide)
            profile_matrix[row][position] += 1

    for position in range(k):
        total = sum(profile_matrix[row][position] for row in range(4))
        for row in range(4):
            profile_matrix[row][position] /= total
    return profile_matrix

def score(motifs, k):
    consensus = ''
    score = 0
    for i in range(k):
        col = [motif[i] for motif in motifs]
        max_freq = max(col, key=col.count)
        consensus += max_freq
        score += sum(1 for base in col if base != max_freq)
    return score, consensus

def highest_prob_kmer(sequence, k, profile_matrix):
    highest_probability_kmer = ''
    highest_prob = -1
    for i in range(len(sequence) - k + 1):
        prob = 1.0
        for j in range(k):
            nucleotide = sequence[i+j]
            if nucleotide == 'A':
                prob *= profile_matrix[0][j]
            elif nucleotide == 'C':
                prob *= profile_matrix[1][j]
            elif nucleotide == 'G':
                prob *= profile_matrix[2][j]
            elif nucleotide == 'T':
                prob *= profile_matrix[3][j]
        if prob > highest_prob:
            highest_prob = prob
            highest_probability_kmer = sequence[i:i+k]
    return highest_probability_kmer

# create motif from profile matrix and sequence 
def motif_from_profile(sequences, k, profile_matrix):
    motifs = []
    for i in range(len(sequences)):
        sequence = sequences[i]
        motifs.append(highest_prob_kmer(sequence, k, profile_matrix))
    return motifs


def randomized_motif(sequences, k, t, N):
    motifs = randomly_select_kmers(sequences, k)
    best_motifs = motifs
    best_score, _ = score(best_motifs, k)
    for _ in range(N):
        # print(_)
        profile = create_profile(motifs, k)
        motifs = motif_from_profile(sequences, k, profile)
        current_score, _ = score(motifs, k)
        if current_score < best_score:
            best_motifs = motifs
            best_score = current_score
    return best_motifs, best_score
   

def run_randomized_motif_scores(filepath, k, N):
    sequences = read_sequences(filepath)
    t = len(sequences)
    best_motifs, best_score = randomized_motif(sequences, k, t, N)
    _, consensus = score(best_motifs, k)
    return best_motifs, best_score, consensus

k = 8
N = 10000
best_motifs, best_score, consensus = run_randomized_motif_scores('hm03.txt', k, N)
print(f"Best motifs: {best_motifs}")
print(f"Score: {best_score}")
print(f"Consensus: {consensus}")