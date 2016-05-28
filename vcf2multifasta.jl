#!/usr/bin/env julia0.3.11
using GZip

### I) functions and libraries
include("merge_vcf.jl")

function correct_ind_names(inds, pop_assign, fnewnames)
	names = ASCIIString[]
	con = open(fnewnames, "w")
	println(con, "CloneName\tNewName")

	for i in 1:length(inds)
		name = "-"
		# dropped individuals get a -
		name = get(pop_assign, inds[i], "-")
		println("found ind $(name)")
		newname = name * @sprintf("_%03d", i)
		push!(names, newname)
		println(con, "$(inds[i])\t$(newname)")
	end
	
	close(con)
	names
end


# merge list of chunks (from ref fasta) and VCF file
function merge(chunks, input, pop_assign, min_dist, config, prefix, fnewnames)
	scaffold = ""
	chunklist = Chunk[]
	first = true
	inds = []

	for line in eachline(input)
        if line != ""
            # comment or header
            if line[1] == '#'
                # header
                if beginswith(line, "#CHROM")
                    # get ind names from header and modify them according to
                    # clone assignment
                    inds = correct_ind_names(split(line)[10:end], pop_assign, fnewnames)
                end
                continue
            end

            # create SNP object
            snp = SNP(strip(line))

            if ! validkind(snp)
                error(coord(snp) * ": invalid type " * snp.kind)
            end

            newchunk = false
            keepmerged = false

            # new scaffold or distance between chunks too big
            # (also triggers on first SNP)
            if snp.scaffold != scaffold || snp.pos > current.chunk.pos_max+min_dist
                scaffold = snp.scaffold
                # we keep chunks as lists per scaffold
                chunklist = chunks[scaffold]
                newchunk = true
            # same scaffold and not too far but still a new chunk
            elseif snp.pos > current.chunk.pos_max
                newchunk = true
                keepmerged = true
                println(coord(snp) * ": continuing with previous chunk")
            end

            println(scaffold, " - ", snp.pos)

            if newchunk
                # we have to find the chunk belonging to this SNP
                p = find_at_pos(chunklist, snp.pos)
                # with a refcall we might still find an overlapping chunk
                if p == 0 && isrefcall(snp)
                    until = endpos(snp)
                    for i in snp.pos+1:until
                        p = find_at_pos(chunklist, p)
                        if p != 0 
                            snp.pos = i
                            break
                        end
                    end
                end
                # if p is still 0 issue a warning, skip this SNP
                if p == 0
                    println(coord(snp) * ": WARNING, position missing in fasta file!")
                    scaffold = ""	# effectively we are between scaffolds here
                    continue
                end

                # log here to avoid spurious messages from missing positions
                if ! keepmerged
                    println("new (merged) target: $(scaffold):$(snp.pos)")
                end
                
                println(coord(snp) * 
                    ": new chunk $p ($(chunklist[p].pos_min)-$(chunklist[p].pos_max))")

                # append new chunk to previous one
                if keepmerged
                    current.chunk = chunklist[p]
                    current.name *= " | " * chunklist[p].id
                # start a new one
                else
                    if !first
                        pad_end(current)
                        save_merged_data(current, prefix, inds)
                    end
                    current = Merged(chunklist[p])
                end
            end

            first = false
            merge(current, snp, config, inds)	# LD: I still don't understand what this line does...
        end

		# end of file
		if eof(input)
			# we won't come here again, so let's save the last one
			pad_end(current)
			save_merged_data(current, prefix, inds)
		end
	end
end

function save_merged_data(merged, prefix, ind_names) 
	# common part of the name
	name_inf = " # $(merged.name) # " * merged.chunk.scaffold * 
		":$(merged.firstpos)-$(merged.chunk.pos_max) # "

	# new file
	fname = prefix * merged.chunk.scaffold * ".$(merged.firstpos).fasta" 
	out = open(fname, "w")

	@assert(!isempty(ind_names) && length(ind_names)*2 == length(merged.data),
		"ERROR wrong number of inds: $(length(ind_names))")

	println("saving $name_inf: $(merged.chunk.pos_max-merged.firstpos+1) -" *
		" $(length(merged.data[1]))")


	order = sortperm(ind_names)

	num = 1
	# go through all haplotypes
	for i in order
		ind_name = isempty(ind_names) ? "" : ind_names[i]

		# these are the unassigned individuals
		if beginswith(ind_name, "-")
			continue
		end

		for str in 0:1
			name = "> $num:$ind_name:$str" * name_inf
			println(out, name)

			num += 1

			hap = merged.data[i*2-1+str]

			# print all sites
			for j in 1:length(hap)
				print(out, hap[j])
				# lines contain 60 letters
				if j%60 == 0
					println(out)
				end
			end
			println(out)
		end
	end

	close(out)
end

function read_clone_assignment(input)
	assign = Dict{ASCIIString, ASCIIString}()

	for line in eachline(input)
		fields = split(line, ',')

		# ignore header line
		if fields[1] == "Clone"
			continue
		end

		# drop marked individuals
		if fields[8] == "F" || fields[8] == "NA"
			continue
		end

		# get rid of _T... suffix
		name = join(split(fields[1], '_')[1:2], "_")
		
		assign[name] = fields[7]
	end

	assign
end


function get_arg(args, i, fun = x->x)
	if length(args) < i+1
		error("expected argument after $(args[i])")
	end

	i += 1
	fun(args[i]), i
end

function print_help()
	println("usage:")
	println("julia multifasta.jlfnewnames <FASTA_REF> <VCF> -a clonefname -N  [[-gqrncdp <ARG>] ...]")  
end


### II) script
if length(ARGS) < 2
	print_help()
	exit()
end

fasta = ARGS[1]
vcf = ARGS[2]

min_dist = 0    
prefix = ""  
clonefname = ""
fnewnames = ""

# config = Config(minGQ, minNR, NACN, minQUAL, minRefQUAL)
config = Config(0, 1, 1.0, 40, 30)

i = 3
while i ≤ length(ARGS)
	arg = ARGS[i]

	if arg == "-g"
		config.minGQ, i = get_arg(ARGS, i, int)
	elseif arg == "-q"
		config.minQUAL, i = get_arg(ARGS, i, int)
	elseif arg == "-r"
		config.minRefQUAL, i = get_arg(ARGS, i, int)
	elseif arg == "-n"
		config.minNR, i = get_arg(ARGS, i, int)
	elseif arg == "-c"
		config.NACN, i = get_arg(ARGS, i, float)
	elseif arg == "-d"
		min_dist, i = get_arg(ARGS, i, int)
	elseif arg == "-p"  
		prefix, i = get_arg(ARGS, i) 
	elseif arg == "-a"
		clonefname, i = get_arg(ARGS, i)
	elseif arg == "-N"	# file to outpout new clone names
		fnewnames, i = get_arg(ARGS, i)
	else
		println("unknown argument $(arg)")
		print_help()
		exit()
	end
	
	i += 1
end

@assert clonefname != "" "please specify a clone assignment file"

pop_assign = open(read_clone_assignment, clonefname)
println(pop_assign)

chunks = open(read_fasta, fasta)

merge(chunks, GZip.open(vcf), pop_assign, min_dist, config, prefix, fnewnames)

