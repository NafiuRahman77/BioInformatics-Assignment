def readinfile():
    # Open the file
    #open an output file
    output = open('out.txt', 'w')
    with open('meme.txt','r') as file:
        
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
print(motifs)