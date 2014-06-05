boundaries = ["periodic"]
filters = ["haar", map(x->string("la",x),[8:2:20]), map(x->string("d",x),[4:2:20]), map(x->string("c",x),[6:6:30])]
levels = [1:4]

include("/Users/tomaskrehlik/Documents/Julia/wavelets/waveletFilters.jl")
include("/Users/tomaskrehlik/Documents/Julia/wavelets/modwt.jl")

data = readcsv("data/testdata.csv")

for i in boundaries, j in filters, l in levels
	run(`Rscript wavelets.R $j $l $i`)
	W = readcsv("W.csv")
	V = readcsv("V.csv")
	(Wjl, Vjl) = modwt(data[:,1], j, l, i)
	@assert all(abs(Wjl - W) .< 2.220446049250313e-12)
	@assert all(abs(Vjl - V) .< 2.220446049250313e-12)
end

run(`rm W.csv`)
run(`rm V.csv`)