export PATH=$HOME/meme/bin:$HOME/meme/libexec/meme-:$PATH
methods=("meme" "streme")
datasets=("hm03" "yst04r" "yst08r")
for m in ${methods[@]}; do
    for d in ${datasets[@]}; do
        for k in {8..15}; do
            mkdir -p "tool_output/$d/$m/k_$k"
            cd "tool_output/$d/$m/k_$k"
            echo -e "\nTool: $m\nDataset: $d\n k=$k\n"
            if [ "$m" = "meme" ]; then
                python3 "../../../../meme_score.py" meme.txt
            elif [ "$m" = "streme" ]; then
                python3 "../../../../streme_score.py" sites.tsv
            fi
            cd ../../../..
        done
    done
done