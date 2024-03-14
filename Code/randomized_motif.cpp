#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include<cstdlib>
#include <chrono>

using namespace std;

double** Profile(vector<string> motifs){
    int total = motifs.size();
    int k = motifs[0].size();
    double **profile = new double*[4];
    for(int i = 0; i < 4; i++){
        profile[i] = new double[k];
    }
    for(int i = 0; i < k; i++){
        int a = 0, c = 0, g = 0, t = 0;
        for(int j = 0; j < total; j++){
            if(motifs[j][i] == 'A'){
                a++;
            } else if(motifs[j][i] == 'C'){
                c++;
            } else if(motifs[j][i] == 'G'){
                g++;
            } else if(motifs[j][i] == 'T'){
                t++;
            }
        }
        profile[0][i] = (double)(a+1) / (total+1); // Laplace's Rule of Succession
        profile[1][i] = (double)(c+1) / (total+1);
        profile[2][i] = (double)(g+1) / (total+1);
        profile[3][i] = (double)(t+1) / (total+1);
    }
    return profile;
}

vector<string> motifs_from_profile(double **profile, vector<string> dna, int k){
    vector<string> motifs;
    int t = dna.size();
    for(int i = 0; i < t; i++){
        string motif;
        double max_prob = 0;
        for(int j = 0; j < dna[0].size() - k; j++){
            string kmer = dna[i].substr(j, k);
            double prob = 1;
            for(int l = 0; l < k; l++){
                if(kmer[l] == 'A'){
                    prob *= profile[0][l];
                } else if(kmer[l] == 'C'){
                    prob *= profile[1][l];
                } else if(kmer[l] == 'G'){
                    prob *= profile[2][l];
                } else if(kmer[l] == 'T'){
                    prob *= profile[3][l];
                }
            }
            if(prob > max_prob){
                max_prob = prob;
                motif = kmer;
            }
        }
        motifs.push_back(motif);
    }
    return motifs;
}

int Score(vector<string> motifs){
    int score = 0;
    int total = motifs.size();
    int k = motifs[0].size();
    for(int i = 0; i < k; i++){
        int a = 0, c = 0, g = 0, t = 0;
        for(int j = 0; j < total; j++){
            if(motifs[j][i] == 'A'){
                a++;
            } else if(motifs[j][i] == 'C'){
                c++;
            } else if(motifs[j][i] == 'G'){
                g++;
            } else if(motifs[j][i] == 'T'){
                t++;
            }
        }
        score += total - max(a, max(c, max(g, t)));
    }
    return score;
}

vector<string> randomized_motif_search(vector<string> dna, int k, int t){
    vector<string> best_motifs;
    vector<string> motifs;
    // Randomly select k-mers Motifs = (Motif1, …, Motift) in each string from Dna
    srand(time(nullptr)); 
    for(int i = 0; i < t; i++){
        int r = rand() % (dna[0].size() - k);
        motifs.push_back(dna[i].substr(r, k));
    }
    // Set BestMotifs to be the motifs
    best_motifs = motifs;
    // while forever
    while(true){
        // Profile ← Profile(Motifs)
        double **profile = Profile(motifs);
        // Motifs ← Motifs(Profile, Dna)
        motifs = motifs_from_profile(profile, dna, k);
        if(Score(motifs) < Score(best_motifs)){
            best_motifs = motifs; 
        } else {
            return best_motifs; 
        }
    }
}

string Consensus(vector<string> motifs){
    string consensus = "";
    int k = motifs[0].size();
    for(int i = 0; i < k; i++){
        int a = 0, c = 0, g = 0, t = 0;
        for(int j = 0; j < motifs.size(); j++){
            if(motifs[j][i] == 'A'){
                a++;
            } else if(motifs[j][i] == 'C'){
                c++;
            } else if(motifs[j][i] == 'G'){
                g++;
            } else if(motifs[j][i] == 'T'){
                t++;
            }
        }
        if(a >= max(c, max(g, t))){
            consensus += 'A';
        } else if(c >= max(a, max(g, t))){
            consensus += 'C';
        } else if(g >= max(a, max(c, t))){
            consensus += 'G';
        } else if(t >= max(a, max(c, g))){
            consensus += 'T';
        }
    }
    return consensus;
}

int main (int argc, char** argv){
    vector<string> dna;
    int t = 0;
    int k = atoi(argv[2]);
    ifstream file(argv[1]);
    string line;
    while(getline(file, line)){
        dna.push_back(line);
        t++;
    }
    file.close();

    auto start = chrono::system_clock::now();
    vector<string> motifs = randomized_motif_search(dna, k, t);
    auto end = chrono::system_clock::now();
    chrono::duration<double> elapsed_seconds = end - start;
    // cout << "Elapsed time: " << elapsed_seconds.count() << "s\n";
    // cout << "Score after 1 iteration: " << Score(motifs) << endl;
    // for(int i = 0; i < motifs.size(); i++){
    //     cout << motifs[i] << endl;
    // }

    vector<string> best_motifs = motifs;
    start = chrono::system_clock::now();
    for(int i = 0; i < 10000; i++){
        motifs = randomized_motif_search(dna, k, t);
        if(Score(motifs) < Score(best_motifs)){
            best_motifs = motifs;
        }
    }
    end = chrono::system_clock::now();
    elapsed_seconds = end - start;
    cout << "Time: " << elapsed_seconds.count() << "\n";
    cout << "Best motifs: ";
    for(int i = 0; i < best_motifs.size(); i++){
        cout << "'" << best_motifs[i] << "'";
        if (i <  best_motifs.size()-1) cout << ", ";
    };
    cout << "]" << endl;
    cout << "Score: " << Score(best_motifs) << endl;
    cout << "Consensus: " << Consensus(best_motifs) << endl;
    return 0;
}