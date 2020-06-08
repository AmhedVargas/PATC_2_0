#!/bin/sh
if [ "$#" -ne 2 ]; then
  echo "Usage: $0 DIRECTORY Threshold" >&2
  exit 1
fi


#Move to temporary directory
cd "${1}" || exit

perl ../../../bin/CreateSeqs2.pl gene.fasta ../../../bin/amino.txt 10000 | ../../../bin/RNAfold | awk -F"-" '{if(NF==1){seq=$0}else{print $2"\t"seq}}' > seqs_energy.txt

cat seqs_energy.txt | perl -pe 's/\(//g' | sort -r -n | head -n "${2}" > sel_seqs.txt
cat sel_seqs.txt| perl -pe 's/U/T/g' |  awk -F"\t" '{print $2}' | perl -pe 's/^AAAA//g' > seqswithoutaaaas.txt



