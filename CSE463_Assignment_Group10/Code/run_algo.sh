# export PATH=$HOME/meme/bin:$HOME/meme/libexec/meme-:$PATH
methods=("gibbs_sampler" "modified_gibbs" "randomized_motif_cpp")
datasets=("hm03" "yst04r" "yst08r")

for d in ${datasets[@]}; do
    for m in ${methods[@]}; do
        for k in {8..15}; do
            echo -e "\nTool: $m\nDataset: $d\nk=$k\n"
            if [ "$m" = "randomized_motif_cpp" ]; then
                ./randomized_motif_cpp "$d.uf.txt" $k
            else
                python3 "$m.py" "$d.uf.txt" $k
            fi
        done
    done
done