# dbcan4 pipeline


# set working directory to project/data
LIST=$(cat 1069MAGS_list.txt)
for i in ${LIST[@]};do
    cd $i
    ln -s ~/THEEND/src/BASH/dbcan4.sh
    cd ..
done

for i in ${LIST[@]};do
    cd $i
    qsub -V dbcan4.sh
    cd ..

done