in=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/encode/encode_macs2_peaks_all/encode_macs2_peaks
out=${in}/beds

mkdir -p $out

for i in ${in}/*.bb; do
	echo $i
	base="$(basename "$i" .bb)"
	echo $base
	bigBedToBed $i ${out}/${base}.bed
done
