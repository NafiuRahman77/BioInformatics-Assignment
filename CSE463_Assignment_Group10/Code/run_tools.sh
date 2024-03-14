export PATH=$HOME/meme/bin:$HOME/meme/libexec/meme-:$PATH
methods=("meme" "streme")
datasets=("hm03" "yst04r" "yst08r")

for d in ${datasets[@]}; do
    for m in ${methods[@]}; do
        for k in {8..15}; do
            mkdir -p "tool_output/$d/$m/k_$k"
            cd "tool_output/$d/$m/k_$k"
            echo -e "\nTool: $m\nDataset: $d\n k=$k\n"
            if [ "$m" = "meme" ]; then
                meme "../../../../$d.txt" -dna -oc . -nostatus -time 14400 -mod zoops -nmotifs 10 -minw $k -maxw $k -objfun classic -revcomp -markov_order 0 
            elif [ "$m" = "streme" ]; then
                streme --verbosity 1 --oc . --dna --totallength 4000000 --time 14400 --minw $k --maxw $k --thresh 0.05 --align center  --p "../../../../$d.txt"
            fi
            cd ../../../..
        done
    done
done