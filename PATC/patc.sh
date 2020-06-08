#!/bin/sh
if [ "$#" -ne 4 ]; then
  echo "Usage: $0 DIRECTORY Threshold Balanced flag-type" >&2
  exit 1
fi


#Move to temporary directory
cd PATC/users/"${1}" || exit

#Initialize output file; cannot be trapped otherwise they will be destroyed upon exit of this shell
echo "" > tab.file
echo "" > DNA-patc.html

#Convert file fasta into a multiline fasta
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' file.fasta | tr "\t" "\n" | fold -w 60 > tmp.file

#Rewrite fasta input
mv tmp.file file.fasta

#Create bat file; first text then parsing to binary
printf "file.fasta\nfile_PATC\n1000000\n%s\n0\n+1\n\"Phasing\"\n\"Date\"\n\"Time\"\n0\n0\n%s\nq\n" "${2}" "${3}"> tmp.file

#Rewwrite bat file
mv tmp.file patc.bat

#Copy compiled PATC; symbolic links does not work cause the Pascal code was compiled having in mind that the bat file was in the same directory
cp ../../NewPATC052115Balanced .

#Run compiled PATC
./NewPATC052115Balanced

#Output only the important info
awk -F"\t" '{OFS="\t"; print $1,$46,$59,$60,$61,$62}' file_PATCpatc.dat | perl -pe 's/\d+:\s//g' > tab.file

if [ "$4" -ne 0 ]; then
	#Modify once again tmp file
	printf "file.fasta\nfile_PATC\n1000000\n%s\n0\n+1\n\"Phasing\"\n\"Date\"\n\"Time\"\n0\n0\nF\nq\n" "${2}" > tmp.file

	#Rewwrite bat file
	mv tmp.file patc.bat

	#Run now with F mode instead balanced
	./NewPATC052115Balanced

	#Output array in fasta
	awk '{if($0 ~ />/){print $0;}else{printf $0}} END{print ""}' file_PATCpatc.num > Filt-array.fasta

	#Old version
	#cat file.fasta | awk '{if($0 ~ />/){print $0;}else{printf $0}} END{print ""}' | cat - Filt-array.fasta | grep -v ">" | awk '{if(NR == 1){i=split($0,array,"")}else{j=split($0,prray,"")}} END{print "<div align=\"justify\"><tt>Potential phasing in DNA input<br><u>N n</u> Nucleotide in phase<br><span style=\"background-color: #FF6347\">A a</span> Adenine belonging to PATC cluster<br><span style=\"background-color: #D4AF37\">T t</span> Thymine belonging to PATC cluster"; for(n=1; n<=i; n++){if(((n-1) % 80) == 0){printf "<br>"};if(prray[n] != "n"){printf "<u>"prray[n]"</u>"}else{printf array[n]}}; print "<tt></div>"}' | perl -pe 's/<u>([Aa])<\/u><u>([Aa])/<u><span style="background-color: #FF6347">$1<\/span>/g' | perl -pe 's/<u>([Tt])<\/u><u>([Tt])/<u><span style="background-color: #D4AF37">$1$2<\/span>/g'  > DNA-patc.html

	#New HTML version
	awk '{if($0 ~ />/){print $0;}else{printf $0}} END{print ""}' file.fasta | cat - Filt-array.fasta | grep -v ">" | awk '{if(NR == 1){i=split($0,array,"")}else{j=split($0,prray,"")}} END{print "<div align=\"justify\"><tt>AT phasing of DNA input<br><u>N n</u> Nucleotide in phase<br><span style=\"background-color: #FF6347\">A a</span> Adenine belonging to PATC cluster<br><span style=\"background-color: #D4AF37\">T t</span> Thymine belonging to PATC cluster"; for(n=1; n<=i; n++){if(((n-1) % 80) == 0){printf "<br>"};if(prray[n] != "n"){printf "<u>"prray[n]"</u>"}else{printf array[n]}}; print "<tt></div>"}' | perl -pe 's/<u>([Aa])/<u><span style="background-color: #FF6347">$1<\/span>/g' | perl -pe 's/<u>([Tt])/<u><span style="background-color: #D4AF37">$1<\/span>/g' > DNA-patc.html
fi

exit
