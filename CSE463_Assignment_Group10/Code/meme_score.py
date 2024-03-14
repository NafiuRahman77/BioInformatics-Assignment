import sys

filename = sys.argv[1]

def score(motifs, k):
    consensus = ''
    score = 0
    for i in range(k):
        col = [motif[i] for motif in motifs]
        max_freq = max(col, key=col.count)
        consensus += max_freq
        score += sum(1 for base in col if base != max_freq)
    return score, consensus

def readinfile():
    # Open the file
    #open an output file
    output = open('out.txt', 'w')
    with open(filename,'r') as file:
        
        start_reading = False
        cnt = 0
        
        #start printing 4 lines after " MEME-1 sites sorted by position p-value"
        for line in file:   
            if "MEME-1 sites sorted by position p-value" in line:
                start_reading = True
                #write to file
                
            if start_reading:
                cnt += 1
                
            if cnt >= 5:
                output.write(line)
                #print(line)
                if "--------------------------------------------------------------------------------" in line:
                    break

    output.close()    
            
def readoutfile():

    motifs=[]

    with open('out.txt','r') as file:
        #cut off the last line

        for line in file:
            if "--------------------------------------------------------------------------------" in line:
                break
            linearr = line.split()
            motifs.append(linearr[5]) 

    return motifs


            
readinfile()

motifs = readoutfile()
print("Best motifs: ", motifs)
scr, cnsns = score(motifs, len(motifs[0]))
print("Score:", scr)
print("Consensus:", cnsns)
