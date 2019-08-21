vcffile="$1"

sed '/^#/ d' $vcffile > "tmp.vcf"
awk '{$NF=""; print $0}' "tmp.vcf" > "tmp2.vcf"
awk '{
for (i=1; i <=NF; i++) {
    if($i ~ /:/) {
        split($i,a,/:/);
        printf "%s\t", a[1];
    }
    else if ($i ~ /;/) {
        printf ".\t";
    }
    else {
       printf "%s\t", $i;
   }
}
printf "\n";
}' tmp2.vcf > $vcffile
rm tmp.vcf tmp2.vcf
