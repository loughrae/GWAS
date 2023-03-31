#conda activate medicc_env - do on command line, not in script
Rscript filter_ascat.R
for aliquot in `cut -f7 ASCAT-forMEDICC.tsv | sort | uniq | grep -v GDC_Aliquot`
do
grep -e $aliquot -e "cn_a" ASCAT-forMEDICC.tsv | cut -f1-6 > ${aliquot}_forMEDICC.tsv
python ../medicc2/medicc2.py ${aliquot}_forMEDICC.tsv ASCAT_MEDICC/ --prefix ${aliquot}_WGD
python ../medicc2/medicc2.py ${aliquot}_forMEDICC.tsv ASCAT_MEDICC/ --prefix ${aliquot}_noWGD --no-wgd
done

