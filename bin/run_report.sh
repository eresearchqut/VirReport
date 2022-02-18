#!/bin/bash

## eResearch Office, QUT
## Created:  31 March 2021
## Last modified: 24 May 2021
## Script: Processes blastN outputs to summarise and report hits to REGULATED and ENDEMIC viruses/viroids.

# Help information to user (i.e., script_name -h or script_name --help)


#Processing tabular files

var=$1
ICTV=$2
echo $var

#STEP0: fetch Top 1 Hits
cat $var | awk '{print $1}' | sort | uniq > ${var}.top1.ids
for i in `cat ${var}.top1.ids`; do echo "fetching top hits..." $i; grep $i $var | head -1 >> ${var}.top1Hits.txt ; done

#STEP1: modify the columns of Galaxy Australia (GA) blast output to the expected format by the BlastTools.jar tool
######  namely: qseqid sgi sacc length pident mismatch gapopen qstart qend qlen sstart send slen sstrand evalue bitscore qcovhsp stitle staxids qseq sseq sseqid qcovs qframe sframe
cat ${var}.top1Hits.txt | sed 's/ /_/g' > ${var}.txt
#STEP2: summarise the GA blastN files
java -jar BlastTools.jar -t blastn ${var}.txt

#filter regulated/edemic/LandPlant
cat summary_${var}.txt | grep "regulated" >> summary_${var}_filtered.txt
cat summary_${var}.txt | grep "endemic" >> summary_${var}_filtered.txt
cat summary_${var}.txt | grep "LandPlant" >> summary_${var}_filtered.txt

#STEP3: fetch unique names from Blast summary reports
cat summary_${var}_filtered.txt | awk '{print $7}' | awk -F "|" '{print $3}'| sort | uniq | sed 's/Species://' > ${var}_uniq.ids

#STEP4: retrieve the best hit for each virus/viroid
echo "processing top hits ..."
for id in `cat ${var}_uniq.ids`
  do
    #fetch the virus name on the summary_blastn file by selecteing longest alignment (column 3) and highest genome coverage (column 5)
    grep $id summary_${var}.txt | sort -k3,3nr -k5,5nr | head -1 >> ${var}_filtered.txt
  done

#print the header of the inital summary_blastn file
cat summary_${var}.txt | head -1 > header

#fetch hits to REGULATED and ENDEMIC viruses
grep "regulated" ${var}_filtered.txt > summary_${var}_REGULATED_viruses_viroids
grep "endemic" ${var}_filtered.txt > summary_${var}_ENDEMIC_viruses_viroids

##### REPORT1 ##### add header to columns
cat header summary_${var}_REGULATED_viruses_viroids > summary_${var}_REGULATED_viruses_viroids.txt
cat header summary_${var}_ENDEMIC_viruses_viroids > summary_${var}_ENDEMIC_viruses_viroids.txt

#fetch genus names of identified hits
awk '{print $7}' summary_${var}_REGULATED_viruses_viroids.txt | awk -F "|" '{print $3}' | sed 's/Species://' | sed 1d > wanted_regulated.names
awk '{print $7}' summary_${var}_ENDEMIC_viruses_viroids.txt | awk -F "|" '{print $3}' | sed 's/Species://' | sed 1d > wanted_endemic.names

#add species to report
paste wanted_regulated.names summary_${var}_REGULATED_viruses_viroids > summary_${var}_REGULATED_viruses_viroids.MOD
paste wanted_endemic.names summary_${var}_ENDEMIC_viruses_viroids > summary_${var}_ENDEMIC_viruses_viroids.MOD

#STEP5: fecth ICTV information
grep -w -F -f wanted_regulated.names $ICTV > wanted_regulated.ICTV
grep -w -F -f wanted_endemic.names $ICTV > wanted_endemic.ICTV

#join reports with ICTV information
join -a 1 -1 1 -2 1 summary_${var}_REGULATED_viruses_viroids.MOD  wanted_regulated.ICTV | tr ' ' '\t' | awk '$4>=70' >  summary_${var}_REGULATED_viruses_viroids_ICTV

#print name of virus/viroid being processed
echo "$id"
join -a 1 -1 1 -2 1 summary_${var}_ENDEMIC_viruses_viroids.MOD  wanted_endemic.ICTV | tr ' ' '\t' | awk '$4>=70' > summary_${var}_ENDEMIC_viruses_viroids_ICTV

#modify header
awk '{print "Species" "\t" $0 "\t" "ICTV_information"}' header > header2

##### REPORT2 ##### add header2 to identified hits
cat header2 summary_${var}_REGULATED_viruses_viroids_ICTV >> summary_${var}_REGULATED_viruses_viroids_ICTV.txt
cat header2 summary_${var}_ENDEMIC_viruses_viroids_ICTV | awk -F"\t" '$1!=""&&$2!=""&&$3!=""' >> summary_${var}_ENDEMIC_viruses_viroids_ICTV.txt

echo "completed!"