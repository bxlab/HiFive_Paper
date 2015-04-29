fends = read.table(fend_fname, header=T)
data = read.table(mat_fname, header=T)

# find bin bounds
minsize = min(fends$coord)
maxsize = max(fends$coord)
binsize = 10000
binstart = floor(minsize / binsize) * binsize
binstop = floor(maxsize / binsize + 1) * binsize
num_bins = (binstop - binstart) / binsize
fend_bins = rep(-1, max(max(data$fend1), max(data$fend2)))
fend_bins[fends$fend] = floor((fends$coord - binstart) / binsize) + 1

# find features
gcc = vector(length=num_bins, mode='numeric')
map = vector(length=num_bins, mode='numeric')
len = vector(length=num_bins, mode='numeric')
size = vector(length=num_bins, mode='numeric')
for(i in 1:nrow(fends)){
	gcc[fend_bins[i]] = gcc[fend_bins[i]] + fends$frag_gc[i]
	map[fend_bins[i]] = map[fend_bins[i]] + fends$map_score[i]
	len[fend_bins[i]] = len[fend_bins[i]] + fends$frag_len[i]
	size[fend_bins[i]] = size[fend_bins[i]] + 1
}
gcc = gcc / max(size, 1)
map = map / max(size, 1)
chr = rep('1', num_bins)
start = (0:(num_bins - 1)) * binsize
end = start + binsize
features = data.frame(start, end, len, gcc, map)

# find data matrix
counts = matrix(nrow=num_bins, ncol=num_bins, data=0)
for(i in 1:nrow(data)){
	bin1 = fend_bins[data$fend1[i]]
	bin2 = fend_bins[data$fend2[i]]
	if(bin1 != bin2 && bin1 != -1 && bin2 != -1){
		counts[bin1, bin2] = counts[bin1, bin2] + data$count[i]
		counts[bin2, bin1] = counts[bin2, bin1] + data$count[i]
	}
}