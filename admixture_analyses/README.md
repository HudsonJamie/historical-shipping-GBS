### Both species follow a similar ADMIXTURE script (the below is for *Microcosmus squamiger*:

- Convert the output .vcf file from ipyrad to .bed:

`vcftools --vcf ./msqua_neutral.vcf --plink --out msqua_neutral`


`plink --file msqua_neutral --make-bed --out ./data/msqua_neutral --noweb`

- Change the location of the `../admixture` based where the admixture program is saved on your computer

`for K in 1 2 3 4 5 6 7 8 9 10; do ../admixture --cv=20 msqua_neutral.bed $K -j40 | tee log_${K}.out; done`

- Extract the cross-validation values for each value of K

`grep -h CV log*.out`

`grep -h CV log*.out > cv_out.csv`

- Zip relevant output files for visualisation using [Clumpak](http://clumpak.tau.ac.il/)  

`zip msqua_neutral.zip sorted_K2_neutral_msqua.Q sorted_K3_neutral_msqua.Q sorted_K4_neutral_msqua.Q`
